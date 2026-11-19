#include "exact_merge_dedup.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <unistd.h>
#include <vector>

namespace {
constexpr size_t kIoBytes = 64 << 10;
constexpr size_t kFanIn = 32;

class TemporaryRuns {
 public:
  explicit TemporaryRuns(const std::string &prefix) {
    std::string pattern = prefix + ".sort-XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back(0);
    if (!mkdtemp(name.data())) throw std::runtime_error("Cannot create sort workspace: " + std::string(std::strerror(errno)));
    path = name.data();
  }
  ~TemporaryRuns() {
    DIR *directory = opendir(path.c_str());
    if (directory) {
      while (auto entry = readdir(directory)) {
        if (std::strcmp(entry->d_name, ".") && std::strcmp(entry->d_name, ".."))
          unlinkat(dirfd(directory), entry->d_name, 0);
      }
      closedir(directory);
    }
    rmdir(path.c_str());
  }
  std::string path;
};

struct Failures {
  std::atomic<bool> failed{false};
  mutable std::mutex mutex;
  std::string message;
  void Save(const std::exception &error) {
    std::lock_guard<std::mutex> lock(mutex);
    if (!failed.load()) { message = error.what(); failed.store(true); }
  }
  void Check() const {
    if (failed.load()) { std::lock_guard<std::mutex> lock(mutex); throw std::runtime_error(message); }
  }
};

struct Counters {
  std::atomic<uint64_t> written{0}, read{0}, runs{0};
};

template <size_t Words>
struct Record {
  std::array<uint64_t, Words> key;
  uint64_t ordinal;
};

template <size_t Words>
class ExternalSorter {
  using Row = Record<Words>;
  struct Slot { std::vector<Row> rows; std::vector<std::string> paths; };
  struct Less {
    bool ordinal;
    bool operator()(const Row &a, const Row &b) const {
      if (ordinal && a.ordinal != b.ordinal) return a.ordinal < b.ordinal;
      if (a.key != b.key) return a.key < b.key;
      return a.ordinal < b.ordinal;
    }
  };
  class Writer {
   public:
    Writer(const std::string &path, Counters &counters) : output_(path, std::ios::binary), counters_(counters) {
      if (!output_) throw std::runtime_error("Cannot open sorted run: " + path);
      rows_.reserve(std::max<size_t>(1, kIoBytes / sizeof(Row)));
      ++counters_.runs;
    }
    void Write(const Row &row) {
      rows_.push_back(row);
      if (rows_.size() == rows_.capacity()) Flush();
    }
    void Close() {
      Flush(); output_.close();
      if (!output_) throw std::runtime_error("Failed to close sorted run");
    }
   private:
    void Flush() {
      if (rows_.empty()) return;
      output_.write(reinterpret_cast<const char *>(rows_.data()), rows_.size() * sizeof(Row));
      if (!output_) throw std::runtime_error("Failed to write sorted run (temporary disk full?)");
      counters_.written.fetch_add(rows_.size() * sizeof(Row));
      rows_.clear();
    }
    std::ofstream output_;
    Counters &counters_;
    std::vector<Row> rows_;
  };
  class Reader {
   public:
    Reader(const std::string &path, Counters &counters, const std::array<uint64_t, Words> *lower = nullptr) : input_(path, std::ios::binary), counters_(counters),
        rows_(std::max<size_t>(1, kIoBytes / sizeof(Row))) {
      if (!input_) throw std::runtime_error("Cannot read sorted run: " + path);
      if (lower) {
        input_.seekg(0, std::ios::end);
        uint64_t first = 0, last = static_cast<uint64_t>(input_.tellg()) / sizeof(Row);
        while (first < last) {
          uint64_t mid = first + (last - first) / 2;
          Row row;
          input_.seekg(mid * sizeof(Row));
          input_.read(reinterpret_cast<char *>(&row), sizeof(Row));
          if (!input_) throw std::runtime_error("Cannot seek sorted run");
          if (row.key < *lower) first = mid + 1; else last = mid;
        }
        input_.seekg(first * sizeof(Row));
      }
    }
    bool Next(Row *row) {
      if (position_ == count_) {
        input_.read(reinterpret_cast<char *>(rows_.data()), rows_.size() * sizeof(Row));
        size_t bytes = static_cast<size_t>(input_.gcount());
        if (input_.bad() || (input_.fail() && !input_.eof()) || bytes % sizeof(Row))
          throw std::runtime_error("Truncated or unreadable sorted run");
        counters_.read.fetch_add(bytes);
        count_ = bytes / sizeof(Row); position_ = 0;
        if (!count_) return false;
      }
      *row = rows_[position_++];
      return true;
    }
   private:
    std::ifstream input_;
    Counters &counters_;
    std::vector<Row> rows_;
    size_t position_ = 0, count_ = 0;
  };

 public:
  ExternalSorter(const std::string &directory, const std::string &name, bool ordinal,
                  uint64_t budget, int workers, Counters &counters, Failures &failures)
      : directory_(directory), name_(name), ordinal_(ordinal), workers_(workers),
        capacity_(std::max<uint64_t>(1, budget / workers / sizeof(Row))),
        counters_(counters), failures_(failures), slots_(workers) {
    static_assert(sizeof(Row) == sizeof(uint64_t) * (Words + 1), "Unexpected sort record padding");
  }

  void Add(const Row &row, int worker) {
    if (failures_.failed.load(std::memory_order_relaxed)) return;
    try {
      auto &rows = slots_[worker].rows;
      if (!rows.capacity()) rows.reserve(capacity_);
      rows.push_back(row);
      if (rows.size() == capacity_) Spill(worker);
    } catch (const std::exception &error) { failures_.Save(error); }
  }

  std::vector<std::string> FinishRuns() {
#pragma omp parallel for schedule(static) num_threads(workers_)
    for (int worker = 0; worker < workers_; ++worker) {
      try { if (!failures_.failed.load()) Spill(worker); }
      catch (const std::exception &error) { failures_.Save(error); }
    }
    failures_.Check();
    std::vector<std::string> paths;
    for (auto &slot : slots_) {
      paths.insert(paths.end(), slot.paths.begin(), slot.paths.end());
      std::vector<Row>().swap(slot.rows);
      slot.paths.clear();
    }
    while (paths.size() > kFanIn) {
      std::vector<std::string> next((paths.size() + kFanIn - 1) / kFanIn);
      const size_t first_id = merge_id_;
      merge_id_ += next.size();
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::min(workers_, 8))
      for (size_t batch = 0; batch < next.size(); ++batch) {
        try {
          size_t begin = batch * kFanIn;
          std::vector<std::string> group(paths.begin() + begin, paths.begin() + std::min(paths.size(), begin + kFanIn));
          std::string path = directory_ + "/" + name_ + "-merge-" + std::to_string(first_id + batch);
          Writer writer(path, counters_);
          Merge(group, [&](const Row &row) { writer.Write(row); });
          writer.Close();
          for (const auto &input : group) Remove(input);
          next[batch] = path;
        } catch (const std::exception &error) { failures_.Save(error); }
      }
      failures_.Check();
      paths.swap(next);
    }
    return paths;
  }

  template <typename Emit>
  void Merge(const std::vector<std::string> &paths, const Emit &emit,
             const std::array<uint64_t, Words> *lower = nullptr,
             const std::array<uint64_t, Words> *upper = nullptr) {
    struct Head { Row row; size_t reader; };
    struct Greater {
      Less less;
      bool operator()(const Head &a, const Head &b) const { return less(b.row, a.row); }
    };
    std::priority_queue<Head, std::vector<Head>, Greater> queue(Greater{Less{ordinal_}});
    std::vector<std::unique_ptr<Reader>> readers;
    for (size_t i = 0; i < paths.size(); ++i) {
      readers.emplace_back(new Reader(paths[i], counters_, lower));
      Row row;
      if (readers.back()->Next(&row)) queue.push({row, i});
    }
    Row previous{};
    bool have_previous = false;
    while (!queue.empty()) {
      Head head = queue.top(); queue.pop();
      if (upper && !(head.row.key < *upper)) break;
      if (ordinal_ || !have_previous || head.row.key != previous.key) {
        emit(head.row); previous = head.row; have_previous = true;
      }
      Row row;
      if (readers[head.reader]->Next(&row)) queue.push({row, head.reader});
    }
  }

  template <typename Emit>
  void MergeParallel(const std::vector<std::string> &paths, const Emit &emit) {
    if (paths.empty()) return;
    // Disjoint key ranges keep equal keys together across all input runs.
    std::vector<std::array<uint64_t, Words>> bounds;
    uint64_t largest = 0;
    std::string sample;
    for (const auto &path : paths) {
      std::ifstream input(path, std::ios::binary | std::ios::ate);
      if (!input) throw std::runtime_error("Cannot sample sorted run");
      uint64_t size = static_cast<uint64_t>(input.tellg());
      if (size > largest) { largest = size; sample = path; }
    }
    const int parallelism = std::min(workers_, 8);
    if (!sample.empty()) {
      std::ifstream input(sample, std::ios::binary);
      for (int i = 1; i < parallelism; ++i) {
        Row row;
        input.seekg((largest / sizeof(Row) * i / parallelism) * sizeof(Row));
        input.read(reinterpret_cast<char *>(&row), sizeof(Row));
        if (!input) throw std::runtime_error("Cannot sample sorted key");
        if (bounds.empty() || bounds.back() != row.key) bounds.push_back(row.key);
      }
    }
#pragma omp parallel for schedule(dynamic, 1) num_threads(parallelism)
    for (size_t part = 0; part <= bounds.size(); ++part) {
      try {
        Merge(paths, [&](const Row &row) { emit(row, static_cast<int>(part)); },
              part ? &bounds[part - 1] : nullptr,
              part < bounds.size() ? &bounds[part] : nullptr);
      } catch (const std::exception &error) { failures_.Save(error); }
    }
    failures_.Check();
  }

  static void Remove(const std::string &path) {
    if (std::remove(path.c_str()) != 0) throw std::runtime_error("Cannot remove consumed sorted run: " + path);
  }

 private:
  void Spill(int worker) {
    auto &slot = slots_[worker];
    if (slot.rows.empty()) return;
    std::sort(slot.rows.begin(), slot.rows.end(), Less{ordinal_});
    std::string path = directory_ + "/" + name_ + "-" + std::to_string(worker) + "-" + std::to_string(slot.paths.size());
    Writer writer(path, counters_);
    const Row *previous = nullptr;
    for (const auto &row : slot.rows) {
      if (ordinal_ || !previous || row.key != previous->key) { writer.Write(row); previous = &row; }
    }
    writer.Close(); slot.paths.push_back(path); slot.rows.clear();
  }
  std::string directory_, name_;
  bool ordinal_;
  int workers_;
  size_t capacity_, merge_id_ = 0;
  Counters &counters_;
  Failures &failures_;
  std::vector<Slot> slots_;
};

template <size_t Words>
class ExactSort : public ExactMergeDedup {
 public:
  ExactSort(int workers, uint64_t budget, const std::string &prefix)
      : directory_(prefix), budget_(budget), workers_(workers), keys_(directory_.path, "keys", false, budget, workers, counters_, failures_) {}
  void AddBase(const uint64_t *key, int worker) override { Add(key, 0, worker); }
  void AddCandidate(const uint64_t *key, uint64_t ordinal) override { Add(key, ordinal, 0); }
  void AddRecord(const uint64_t *key, uint64_t priority, int worker) override { Add(key, priority, worker); }
  void FinishUnordered(const std::function<void(const uint64_t *, uint64_t, int)> &emit) override {
    auto paths = keys_.FinishRuns();
    keys_.MergeParallel(paths, [&](const Record<Words> &row, int worker) {
      emit(row.key.data(), row.ordinal, worker);
    });
    for (const auto &path : paths) ExternalSorter<Words>::Remove(path);
  }
  void Check() const override { failures_.Check(); }
  void Finish(const std::function<void(const uint64_t *, uint64_t)> &emit) override {
    auto paths = keys_.FinishRuns();
    ExternalSorter<Words> order(directory_.path, "order", true, budget_, workers_, counters_, failures_);
    // A base record has ordinal zero, so it wins every equal-key group. All
    // other groups retain the first generated candidate, then regain that order.
    keys_.MergeParallel(paths, [&](const Record<Words> &row, int worker) { if (row.ordinal) order.Add(row, worker); });
    failures_.Check();
    for (const auto &path : paths) ExternalSorter<Words>::Remove(path);
    auto accepted = order.FinishRuns();
    order.Merge(accepted, [&](const Record<Words> &row) { emit(row.key.data(), row.ordinal); });
    for (const auto &path : accepted) ExternalSorter<Words>::Remove(path);
  }
  ExactDedupStats Stats() const override {
    return {counters_.written.load(), counters_.read.load(), counters_.runs.load(), budget_};
  }
 private:
  void Add(const uint64_t *key, uint64_t ordinal, int worker) {
    Record<Words> row;
    std::copy_n(key, Words, row.key.begin()); row.ordinal = ordinal;
    keys_.Add(row, worker);
  }
  TemporaryRuns directory_;
  uint64_t budget_;
  int workers_;
  Counters counters_;
  Failures failures_;
  ExternalSorter<Words> keys_;
};
}  // namespace

std::unique_ptr<ExactMergeDedup> MakeExactMergeDedup(unsigned words, int workers,
                                                  uint64_t buffer_bytes,
                                                  const std::string &prefix) {
  if (workers < 1 || buffer_bytes < (64ULL << 10))
    throw std::runtime_error("Invalid exact-sort worker count or record-buffer budget");
#define SORT_CASE(N) case N: return std::unique_ptr<ExactMergeDedup>(new ExactSort<N>(workers, buffer_bytes, prefix))
  switch (words) {
    SORT_CASE(1); SORT_CASE(2); SORT_CASE(3); SORT_CASE(4);
    SORT_CASE(5); SORT_CASE(6); SORT_CASE(7); SORT_CASE(8);
    default: throw std::runtime_error("Unsupported packed width for exact merge sort");
  }
#undef SORT_CASE
}
