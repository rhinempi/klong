#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "definitions.h"
#include "kmlib/kmcompactvector.h"
#include "utils/utils.h"
#include "sequence/io/edge/edge_io_meta.h"
#include "utils/options_description.h"
#include "end_fm_index.h"
#include "exact_merge_dedup.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr uint64_t kDefaultMaxBloomBytes = 60ULL * 1024 * 1024 * 1024;
constexpr uint64_t kSplitMixSeed = 0x9e3779b97f4a7c15ULL;
constexpr size_t kEdgeBufferBytes = 1ULL << 20;

struct Options {
  std::string mode{"reduce"};
  std::string short_prefix;
  std::string long_prefix;
  std::string output_prefix;
  std::string contig_file;
  std::string long_contig_file;
  unsigned short_k{0};
  unsigned long_k{0};
  int num_threads{1};
  double max_bloom_gb{60.0};
  double reserve_gb{1.0};
  int max_hashes{16};
  double host_mem{0};
  double reduce_fpr{1e-15};
  int human_readable{0};
  std::string merge_index{"hash"};
  std::string merge_dedup{"sort"};
  unsigned contig_edge_multiplicity{50};
  bool reduce_dedup{false};
  unsigned read_counted_shards{0};
  double sort_memory_mb{128};
};

struct KmerKey {
  static constexpr int kWords = (kMaxK + 1 + 31) / 32;
  uint64_t words[kWords];

  KmerKey() { std::fill(words, words + kWords, 0); }

  bool operator==(const KmerKey &other) const {
    return std::memcmp(words, other.words, sizeof(words)) == 0;
  }
};

uint64_t SplitMix64(uint64_t x) {
  x += kSplitMixSeed;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t HashKey(const KmerKey &key, uint64_t seed, int words = KmerKey::kWords) {
  uint64_t h = seed;
  for (int i = 0; i < words; ++i) {
    h = SplitMix64(h ^ key.words[i]);
  }
  return h;
}

struct KmerKeyHash {
  size_t operator()(const KmerKey &key) const {
    return static_cast<size_t>(HashKey(key, 0x9e3779b97f4a7c15ULL));
  }
};

struct LongEdgeOrigin {
  std::string edge;
  std::string side;
};

unsigned EdgeLength(const EdgeIoMetadata &meta) { return meta.kmer_size + 1; }

uint8_t GetBase(const uint32_t *edge, unsigned index) {
  unsigned word = index / kCharsPerEdgeWord;
  unsigned offset = index % kCharsPerEdgeWord;
  unsigned shift = (kCharsPerEdgeWord - 1 - offset) * kBitsPerEdgeChar;
  return static_cast<uint8_t>((edge[word] >> shift) & kEdgeCharMask);
}

void SetPackedBase(KmerKey *key, unsigned index, uint8_t base) {
  unsigned bit = index * 2;
  unsigned word = bit / 64;
  unsigned shift = 62 - (bit % 64);
  key->words[word] |= static_cast<uint64_t>(base) << shift;
}

uint8_t PackedKeyBase(const KmerKey &key, unsigned index) {
  return (key.words[index / 32] >> (62 - 2 * (index % 32))) & 3;
}

bool PackedKeyLess(const KmerKey &lhs, const KmerKey &rhs) {
  for (int i = 0; i < KmerKey::kWords; ++i) {
    if (lhs.words[i] != rhs.words[i]) return lhs.words[i] < rhs.words[i];
  }
  return false;
}

KmerKey PackedKeySlice(const KmerKey &key, unsigned offset, unsigned length) {
  KmerKey result;
  unsigned shift = 2 * (offset % 32);
  unsigned first_word = offset / 32;
  unsigned words = (length + 31) / 32;
  for (unsigned i = 0; i < words; ++i) {
    result.words[i] = key.words[first_word + i] << shift;
    if (shift && first_word + i + 1 < KmerKey::kWords) {
      result.words[i] |= key.words[first_word + i + 1] >> (64 - shift);
    }
  }
  if (length % 32) result.words[words - 1] &= ~uint64_t{0} << (64 - 2 * (length % 32));
  return result;
}

uint64_t ReverseComplementWord(uint64_t word) {
  uint64_t x = ~word;
  x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
  x = ((x >> 4) & 0x0f0f0f0f0f0f0f0fULL) | ((x & 0x0f0f0f0f0f0f0f0fULL) << 4);
  x = ((x >> 8) & 0x00ff00ff00ff00ffULL) | ((x & 0x00ff00ff00ff00ffULL) << 8);
  x = ((x >> 16) & 0x0000ffff0000ffffULL) | ((x & 0x0000ffff0000ffffULL) << 16);
  return (x >> 32) | (x << 32);
}

KmerKey PackedReverseComplement(const KmerKey &key, unsigned length) {
  KmerKey reversed;
  unsigned words = (length + 31) / 32;
  for (unsigned i = 0; i < words; ++i)
    reversed.words[i] = ReverseComplementWord(key.words[words - 1 - i]);
  return PackedKeySlice(reversed, words * 32 - length, length);
}

class MergeQueryIndex {
 public:
  MergeQueryIndex(unsigned length, const std::string &backend)
      : compact_(length <= 32), use_fm_(backend == "fm"), fm_(length) {}

  void Register(const KmerKey &key) {
    if (use_fm_) { fm_.Add(key.words); return; }
    bool added = compact_ ? small_.emplace(key.words[0], hits_.size()).second
                          : large_.emplace(key, hits_.size()).second;
    if (added) hits_.emplace_back();
  }

  void Build(int threads) {
    if (use_fm_) { fm_.Build(threads); hits_.resize(fm_.size()); }
  }

  size_t FindId(const KmerKey &key) const {
    if (use_fm_) return fm_.Find(key.words);
    if (compact_) {
      auto it = small_.find(key.words[0]);
      return it == small_.end() ? EndFmIndex::Missing() : it->second;
    }
    auto it = large_.find(key);
    return it == large_.end() ? EndFmIndex::Missing() : it->second;
  }

  std::vector<KmerKey> *Find(const KmerKey &key) {
    size_t id = FindId(key);
    return id == EndFmIndex::Missing() ? nullptr : &hits_[id];
  }

  std::vector<KmerKey> &Hits(size_t id) { return hits_[id]; }
  size_t size() const { return hits_.size(); }
  size_t fm_bytes() const { return use_fm_ ? fm_.bytes() : 0; }

 private:
  bool compact_, use_fm_;
  EndFmIndex fm_;
  std::vector<std::vector<KmerKey>> hits_;
  std::unordered_map<uint64_t, size_t> small_;
  std::unordered_map<KmerKey, size_t, KmerKeyHash> large_;
};

struct PackedMergeContig {
  kmlib::CompactVector<2, uint64_t, kmlib::kBigEndian> bases;
  std::vector<size_t> invalid;
  std::string name;

  PackedMergeContig(const std::string &sequence, const std::string &header,
                    bool diagnostic) : name(diagnostic ? header + "/rev" : "") {
    bases.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
      uint8_t base = 0;
      switch (*it) {
        case 'A': case 'a': break;
        case 'C': case 'c': base = 1; break;
        case 'G': case 'g': base = 2; break;
        case 'T': case 't': base = 3; break;
        default: invalid.push_back(bases.size()); break;
      }
      bases.push_back(base);
    }
  }

  uint8_t Base(size_t offset) const {
    if (!invalid.empty() && std::binary_search(invalid.begin(), invalid.end(), offset)) return 4;
    return bases[offset];
  }

  bool EndKey(size_t offset, unsigned length, KmerKey *key) const {
    *key = KmerKey{};
    for (unsigned i = 0; i < length; ++i) {
      uint8_t base = Base(offset + i);
      if (base > 3) return false;
      SetPackedBase(key, i, base);
    }
    return true;
  }
};

struct PackedMergeView {
  const PackedMergeContig &contig;
  const KmerKey *edge;
  size_t split;
  size_t tail_offset;
  bool edge_first;
  size_t length;

  uint8_t operator()(size_t offset) const {
    if (edge == nullptr) return contig.Base(offset);
    if (edge_first) {
      return offset < split ? PackedKeyBase(*edge, offset)
                            : contig.Base(tail_offset + offset - split);
    }
    return offset < split ? contig.Base(offset) : PackedKeyBase(*edge, offset - split);
  }
};

template <typename BaseAt>
bool CanonicalPackedWindow(const BaseAt &base_at, size_t offset, unsigned length,
                            KmerKey *result) {
  KmerKey forward, reverse;
  for (unsigned i = 0; i < length; ++i) {
    uint8_t base = base_at(offset + i);
    if (base > 3) return false;
    SetPackedBase(&forward, i, base);
    SetPackedBase(&reverse, length - 1 - i, 3 - base);
  }
  *result = PackedKeyLess(reverse, forward) ? reverse : forward;
  return true;
}

int CompareForwardAndReverseComplement(const uint32_t *edge, unsigned offset,
                                       unsigned len) {
  for (unsigned i = 0; i < len; ++i) {
    uint8_t f = GetBase(edge, offset + i);
    uint8_t r = 3 - GetBase(edge, offset + len - 1 - i);
    if (f != r) return f < r ? -1 : 1;
  }
  return 0;
}

KmerKey CanonicalKeyFromEdge(const uint32_t *edge, unsigned offset,
                             unsigned len) {
  KmerKey key;
  bool use_forward = CompareForwardAndReverseComplement(edge, offset, len) <= 0;
  for (unsigned i = 0; i < len; ++i) {
    uint8_t base = use_forward ? GetBase(edge, offset + i)
                               : static_cast<uint8_t>(3 - GetBase(edge, offset + len - 1 - i));
    SetPackedBase(&key, i, base);
  }
  return key;
}

KmerKey RawKeyFromEdge(const uint32_t *edge, unsigned offset, unsigned len) {
  KmerKey key;
  // Repack 32-bit edge words into 64-bit windows without visiting each base.
  for (unsigned i = 0; i < (len + 31) / 32; ++i) {
    unsigned start = offset + i * 32;
    unsigned length = std::min(32U, len - i * 32);
    unsigned word = start / 16, skip = start % 16, shift = skip * 2;
    uint64_t packed = uint64_t{edge[word]} << 32;
    if (skip + length > 16) packed |= edge[word + 1];
    packed <<= shift;
    if (skip + length > 32) packed |= uint64_t{edge[word + 2]} >> (32 - shift);
    if (length < 32) packed &= ~uint64_t{0} << (64 - 2 * length);
    key.words[i] = packed;
  }
  return key;
}

KmerKey RawReverseComplementKeyFromEdge(const uint32_t *edge, unsigned offset,
                                        unsigned len) {
  return PackedReverseComplement(RawKeyFromEdge(edge, offset, len), len);
}

std::string ReverseComplement(const std::string &seq) {
  std::string rc;
  rc.reserve(seq.size());
  for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
    switch (*it) {
      case 'A': case 'a': rc.push_back('T'); break;
      case 'C': case 'c': rc.push_back('G'); break;
      case 'G': case 'g': rc.push_back('C'); break;
      case 'T': case 't': rc.push_back('A'); break;
      default: rc.push_back('N'); break;
    }
  }
  return rc;
}

std::string BasesFromEdge(const uint32_t *edge, unsigned offset, unsigned len) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  std::string s;
  s.reserve(len);
  for (unsigned i = 0; i < len; ++i) s.push_back(bases[GetBase(edge, offset + i)]);
  return s;
}

std::string BasesFromKey(const KmerKey &key, unsigned len) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  std::string s;
  s.reserve(len);
  for (unsigned i = 0; i < len; ++i) {
    unsigned bit = i * 2;
    unsigned word = bit / 64;
    unsigned shift = 62 - (bit % 64);
    s.push_back(bases[(key.words[word] >> shift) & 3]);
  }
  return s;
}

uint32_t Multiplicity(const uint32_t *edge, uint32_t words_per_edge) {
  return edge[words_per_edge - 1] & kMaxMul;
}

void PackKeyAsEdge(const KmerKey &key, unsigned len, uint32_t words_per_edge,
                   uint32_t multiplicity, std::vector<uint32_t> *edge) {
  edge->assign(words_per_edge, 0);
  for (unsigned i = 0; i < len; ++i) {
    unsigned bit64 = i * 2;
    unsigned word64 = bit64 / 64;
    unsigned shift64 = 62 - (bit64 % 64);
    uint32_t base = (key.words[word64] >> shift64) & 3;
    unsigned word32 = i / kCharsPerEdgeWord;
    unsigned offset32 = i % kCharsPerEdgeWord;
    unsigned shift32 = (kCharsPerEdgeWord - 1 - offset32) * kBitsPerEdgeChar;
    (*edge)[word32] |= base << shift32;
  }
  (*edge)[words_per_edge - 1] |= std::min<uint32_t>(multiplicity, kMaxMul);
}

EdgeIoMetadata LoadMetadata(const std::string &prefix) {
  std::ifstream is(prefix + ".edges.info");
  if (!is) xfatal("Cannot open {s}.edges.info\n", prefix.c_str());
  EdgeIoMetadata meta;
  meta.Deserialize(is);
  return meta;
}

uint64_t FileRecordCount(const std::string &prefix, const EdgeIoMetadata &meta,
                         unsigned file_id) {
  if (!meta.is_sorted) {
    std::ifstream is(prefix + ".edges." + std::to_string(file_id),
                     std::ifstream::binary | std::ifstream::ate);
    if (!is) xfatal("Cannot open {s}.edges.{}\n", prefix.c_str(), file_id);
    uint64_t bytes = static_cast<uint64_t>(is.tellg());
    uint64_t stride = sizeof(uint32_t) * meta.words_per_edge;
    if (stride == 0 || bytes % stride) xfatal("Invalid edge shard size for {s}.edges.{}\n", prefix.c_str(), file_id);
    return bytes / stride;
  }
  uint64_t n = 0;
  for (const auto &bucket : meta.buckets) {
    if (bucket.file_id == static_cast<int>(file_id)) {
      n += static_cast<uint64_t>(bucket.total_number);
    }
  }
  return n;
}

uint64_t TotalRecords(const EdgeIoMetadata &meta) {
  if (!meta.is_sorted) return static_cast<uint64_t>(meta.num_edges);
  uint64_t n = 0;
  for (const auto &bucket : meta.buckets) n += static_cast<uint64_t>(bucket.total_number);
  return n;
}

uint64_t AvailableMemoryBytes() {
  std::ifstream is("/proc/meminfo");
  std::string key;
  uint64_t value_kb;
  std::string unit;
  uint64_t available = 0;
  while (is >> key >> value_kb >> unit) {
    if (key == "MemAvailable:") { available = value_kb * 1024ULL; break; }
  }
  // MemAvailable describes the host, not necessarily the process's container.
  auto apply_limit = [&](const std::string &limit_path, const std::string &usage_path) {
    uint64_t limit, usage;
    std::ifstream limit_file(limit_path), usage_file(usage_path);
    if (limit_file >> limit && usage_file >> usage) {
      uint64_t remaining = limit > usage ? limit - usage : 0;
      if (!available || remaining < available) available = remaining;
      if (remaining == 0) available = 1;
    }
  };
  apply_limit("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current");
  apply_limit("/sys/fs/cgroup/memory/memory.limit_in_bytes",
              "/sys/fs/cgroup/memory/memory.usage_in_bytes");
  std::ifstream groups("/proc/self/cgroup");
  std::string group;
  while (std::getline(groups, group)) {
    auto first = group.find(':'), last = group.rfind(':');
    if (first == std::string::npos || last == first) continue;
    std::string path = group.substr(last + 1);
    if (path.empty() || path[0] != '/' || path.find("..") != std::string::npos) continue;
    if (group.substr(0, last) == "0:") {
      apply_limit("/sys/fs/cgroup" + path + "/memory.max",
                  "/sys/fs/cgroup" + path + "/memory.current");
    } else {
      std::string controllers = "," + group.substr(first + 1, last - first - 1) + ",";
      if (controllers.find(",memory,") != std::string::npos)
        apply_limit("/sys/fs/cgroup/memory" + path + "/memory.limit_in_bytes",
                    "/sys/fs/cgroup/memory" + path + "/memory.usage_in_bytes");
    }
  }
  return available;
}

class BloomFilter {
 public:
  BloomFilter(uint64_t num_bits, int num_hashes,
              int key_words = KmerKey::kWords, int threads = 1)
      : num_bits_(std::max<uint64_t>(64, num_bits)),
        num_words_((num_bits_ + 63) / 64),
        num_hashes_(std::max(1, num_hashes)),
        key_words_(key_words),
        bit_mask_((num_bits_ & (num_bits_ - 1)) == 0 ? num_bits_ - 1 : 0),
        words_(new std::atomic<uint64_t>[num_words_]) {
#pragma omp parallel for schedule(static) num_threads(threads) if(num_words_ >= 131072)
    for (uint64_t i = 0; i < num_words_; ++i) words_[i].store(0, std::memory_order_relaxed);
  }

  void Add(const KmerKey &key) {
    uint64_t h1 = HashKey(key, 0x243f6a8885a308d3ULL, key_words_);
    uint64_t h2 = HashKey(key, 0x13198a2e03707344ULL, key_words_) | 1ULL;
    for (int i = 0; i < num_hashes_; ++i) {
      uint64_t hash = h1 + static_cast<uint64_t>(i) * h2;
      uint64_t bit = bit_mask_ ? hash & bit_mask_ : hash % num_bits_;
      words_[bit >> 6].fetch_or(1ULL << (bit & 63), std::memory_order_relaxed);
    }
  }

  bool Contains(const KmerKey &key) const {
    uint64_t h1 = HashKey(key, 0x243f6a8885a308d3ULL, key_words_);
    uint64_t h2 = HashKey(key, 0x13198a2e03707344ULL, key_words_) | 1ULL;
    for (int i = 0; i < num_hashes_; ++i) {
      uint64_t hash = h1 + static_cast<uint64_t>(i) * h2;
      uint64_t bit = bit_mask_ ? hash & bit_mask_ : hash % num_bits_;
      if ((words_[bit >> 6].load(std::memory_order_relaxed) & (1ULL << (bit & 63))) == 0) {
        return false;
      }
    }
    return true;
  }

  uint64_t num_bits() const { return num_bits_; }
  uint64_t num_bytes() const { return num_words_ * sizeof(uint64_t); }
  int num_hashes() const { return num_hashes_; }

 private:
  uint64_t num_bits_;
  uint64_t num_words_;
  int num_hashes_;
  int key_words_;
  uint64_t bit_mask_;
  std::unique_ptr<std::atomic<uint64_t>[]> words_;
};

class EdgeOutputWriter {
 public:
  EdgeOutputWriter(const std::string &prefix, uint32_t k, uint32_t words_per_edge,
                   int num_files, bool human_readable)
      : prefix_(prefix), k_(k), words_per_edge_(words_per_edge),
        human_readable_(human_readable), counts_(std::max(1, num_files), 0),
        buffers_(std::max(1, num_files)) {
    int n = std::max(1, num_files);
    for (int i = 0; i < n; ++i) {
      binary_.emplace_back(new std::ofstream(prefix + ".edges." + std::to_string(i),
                                             std::ofstream::binary | std::ofstream::out));
      if (!*binary_.back()) xfatal("Cannot open output edge shard {}\n", i);
      if (human_readable_) {
        text_.emplace_back(new std::ofstream(prefix + ".edges." + std::to_string(i) + ".txt"));
        *text_.back() << "#edge\tmultiplicity\tsource\n";
      }
    }
  }

  void Write(const uint32_t *edge, int file_id, const char *source) {
    file_id %= static_cast<int>(binary_.size());
    auto &buffer = buffers_[file_id];
    if (buffer.capacity() == 0) buffer.reserve(kEdgeBufferBytes / sizeof(uint32_t) + words_per_edge_);
    buffer.insert(buffer.end(), edge, edge + words_per_edge_);
    if (buffer.size() * sizeof(uint32_t) >= kEdgeBufferBytes) Flush(file_id);
    if (human_readable_) {
      *text_[file_id] << BasesFromEdge(edge, 0, k_ + 1) << '\t'
                      << Multiplicity(edge, words_per_edge_) << '\t' << source << '\n';
    }
    ++counts_[file_id];
  }

  void Finalize() {
    for (size_t i = 0; i < binary_.size(); ++i) {
      Flush(i);
      binary_[i]->close();
      if (!*binary_[i]) xfatal("Failed to close output edge shard {}\n", i);
    }
    for (auto &f : text_) f->close();
    EdgeIoMetadata meta;
    meta.kmer_size = k_;
    meta.words_per_edge = words_per_edge_;
    meta.num_files = counts_.size();
    meta.num_edges = 0;
    meta.is_sorted = false;
    for (auto c : counts_) meta.num_edges += c;
    std::ofstream info(prefix_ + ".edges.info");
    meta.Serialize(info);
    info.close();
    if (!info) xfatal("Failed to write edge metadata for {s}\n", prefix_.c_str());
  }

 private:
  void Flush(size_t file_id) {
    auto &buffer = buffers_[file_id];
    if (buffer.empty()) return;
    binary_[file_id]->write(reinterpret_cast<const char *>(buffer.data()),
                           buffer.size() * sizeof(uint32_t));
    if (!*binary_[file_id]) xfatal("Failed to write output edge shard {}\n", file_id);
    buffer.clear();
  }
  std::string prefix_;
  uint32_t k_;
  uint32_t words_per_edge_;
  bool human_readable_;
  std::vector<std::unique_ptr<std::ofstream>> binary_;
  std::vector<std::unique_ptr<std::ofstream>> text_;
  std::vector<int64_t> counts_;
  std::vector<std::vector<uint32_t>> buffers_;
};

template <typename Fn>
void ForEachEdgeRange(const std::string &prefix, const EdgeIoMetadata &meta,
                     unsigned file_id, uint64_t first, uint64_t n, const Fn &fn) {
  if (n == 0) return;

  std::ifstream is(prefix + ".edges." + std::to_string(file_id),
                   std::ifstream::binary | std::ifstream::in);
  if (!is) xfatal("Cannot open {s}.edges.{}\n", prefix.c_str(), file_id);

  if (meta.words_per_edge == 0) xfatal("Invalid zero-sized edge record\n");
  is.seekg(first * sizeof(uint32_t) * meta.words_per_edge);
  uint64_t block_records = std::max<uint64_t>(1, kEdgeBufferBytes /
                                               (sizeof(uint32_t) * meta.words_per_edge));
  std::vector<uint32_t> block(std::min(n, block_records) * meta.words_per_edge);
  for (uint64_t begin = 0; begin < n; begin += block_records) {
    uint64_t count = std::min(block_records, n - begin);
    is.read(reinterpret_cast<char *>(block.data()),
            count * sizeof(uint32_t) * meta.words_per_edge);
    if (!is) xfatal("Unexpected EOF reading {s}.edges.{}\n", prefix.c_str(), file_id);
    for (uint64_t i = 0; i < count; ++i) fn(block.data() + i * meta.words_per_edge);
  }
}

template <typename Fn>
void ForEachEdgeInFile(const std::string &prefix, const EdgeIoMetadata &meta,
                       unsigned file_id, const Fn &fn) {
  ForEachEdgeRange(prefix, meta, file_id, 0, FileRecordCount(prefix, meta, file_id), fn);
}

struct EdgeRange {
  unsigned file_id;
  uint64_t first;
  uint64_t count;
};

std::vector<EdgeRange> ReductionRanges(const std::string &prefix,
                                      const EdgeIoMetadata &meta, int threads) {
  if (meta.words_per_edge == 0) xfatal("Invalid zero-sized edge record\n");
  uint64_t stride = sizeof(uint32_t) * meta.words_per_edge;
  uint64_t tasks = static_cast<uint64_t>(threads) * 4;
  uint64_t records = TotalRecords(meta);
  uint64_t chunk = std::max<uint64_t>(std::max<uint64_t>(1, kEdgeBufferBytes / stride),
                                      records / tasks + (records % tasks != 0));
  std::vector<EdgeRange> ranges;
  for (unsigned file = 0; file < meta.num_files; ++file) {
    uint64_t count = FileRecordCount(prefix, meta, file);
    for (uint64_t first = 0; first < count; first += chunk)
      ranges.push_back({file, first, std::min(chunk, count - first)});
  }
  return ranges;
}

int ChooseHashCount(uint64_t bits, uint64_t items, int max_hashes) {
  if (items == 0) return 1;
  double optimal = (static_cast<double>(bits) / static_cast<double>(items)) * std::log(2.0);
  return static_cast<int>(std::max(1.0, std::min<double>(max_hashes, std::round(optimal))));
}

double EstimateFpr(uint64_t bits, uint64_t items, int hashes) {
  if (bits == 0) return 1.0;
  double exponent = -static_cast<double>(hashes) * static_cast<double>(items) /
                    static_cast<double>(bits);
  return std::pow(-std::expm1(exponent), hashes);
}

uint64_t SelectBloomByteCap(const Options &opt, uint64_t io_bytes = 0) {
  uint64_t requested = static_cast<uint64_t>(opt.max_bloom_gb * 1024.0 * 1024.0 * 1024.0);
  uint64_t reserve = static_cast<uint64_t>(opt.reserve_gb * 1024.0 * 1024.0 * 1024.0);
  uint64_t available = AvailableMemoryBytes();
  if (opt.max_bloom_gb == 0) requested = kDefaultMaxBloomBytes;
  if (requested < 64) xfatal("Bloom memory cap must allow at least 64 bytes\n");
  if (opt.host_mem) available = available ? std::min(available, static_cast<uint64_t>(opt.host_mem))
                                        : static_cast<uint64_t>(opt.host_mem);
  if (available == 0) return requested;
  reserve = std::min(reserve, available / 8);
  if (available <= reserve || available - reserve <= io_bytes + 64)
    xfatal("Memory budget is too small for reduction I/O buffers; reduce thread count or raise --host_mem\n");
  return std::min(requested, available - reserve - io_bytes);
}

struct ReductionBloomShape {
  uint64_t bytes;
  int hashes;
  double estimated_fpr;
};

ReductionBloomShape SelectReductionBloom(uint64_t items, uint64_t cap, const Options &opt) {
  if (cap < 64) xfatal("Bloom memory cap must allow at least 64 bytes\n");
  items = std::max<uint64_t>(1, items);
  int preferred = std::min(12, opt.max_hashes);
  double bits = -static_cast<double>(preferred) * items /
                std::log1p(-std::pow(opt.reduce_fpr, 1.0 / preferred));
  uint64_t desired = static_cast<uint64_t>(std::min<double>(cap, std::ceil(bits / 8)));
  desired = std::min(cap, std::max<uint64_t>(1ULL << 20, desired));
  uint64_t bytes = 8;
  while (bytes <= cap / 2 && bytes < desired) bytes *= 2;
  if (bytes < desired) bytes = cap & ~uint64_t{7};
  // A power-of-two bit array replaces integer division with a mask per probe.
  int optimal = ChooseHashCount(bytes * 8, items, opt.max_hashes);
  int hashes = optimal;
  for (int h = 1; h <= optimal; ++h) {
    if (EstimateFpr(bytes * 8, items, h) <= opt.reduce_fpr) { hashes = h; break; }
  }
  return {bytes, hashes, EstimateFpr(bytes * 8, items, hashes)};
}

bool ReadNextFasta(std::ifstream *in, std::string *name, std::string *seq) {
  static std::string pending_name;
  seq->clear();
  std::string line;
  if (!pending_name.empty()) {
    *name = pending_name;
    pending_name.clear();
  } else {
    name->clear();
    while (std::getline(*in, line)) {
      if (!line.empty() && line[0] == '>') {
        *name = line.substr(1);
        break;
      }
    }
  }
  if (name->empty()) return false;
  while (std::getline(*in, line)) {
    if (!line.empty() && line[0] == '>') {
      pending_name = line.substr(1);
      return true;
    }
    seq->append(line);
  }
  return !seq->empty();
}

Options ParseOptions(int argc, char **argv) {
  Options opt;
  OptionsDescription desc;
  desc.AddOption("mode", "M", opt.mode, "reduce, merge, or index-bench");
  desc.AddOption("short_prefix", "s", opt.short_prefix, "prefix of shorter k binary edge files");
  desc.AddOption("long_prefix", "l", opt.long_prefix, "prefix of longer/current k binary edge files");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix for edge files");
  desc.AddOption("contig", "c", opt.contig_file, "short/pruned contigs used by merge mode overlap extension");
  desc.AddOption("long_contig", "C", opt.long_contig_file, "long/pruned contigs scanned directly for long-k edges in merge mode");
  desc.AddOption("short_k", "k", opt.short_k, "shorter k-mer size");
  desc.AddOption("long_k", "K", opt.long_k, "longer/current k-mer size");
  desc.AddOption("num_cpu_threads", "t", opt.num_threads, "number of CPU threads/output files");
  desc.AddOption("max_bloom_gb", "m", opt.max_bloom_gb, "max Bloom filter memory in GB");
  desc.AddOption("reserve_gb", "r", opt.reserve_gb, "RAM to leave unused in GB");
  desc.AddOption("max_hashes", "H", opt.max_hashes, "maximum number of Bloom hash probes");
  desc.AddOption("host_mem", "", opt.host_mem, "total memory budget in bytes (0: available memory)");
  desc.AddOption("reduce_fpr", "", opt.reduce_fpr, "target estimated false-positive rate for reduction");
  desc.AddOption("human_readable", "", opt.human_readable, "write human-readable debug reports");
  desc.AddOption("merge_index", "", opt.merge_index, "exact contig-end index: fm or hash");
  desc.AddOption("merge_dedup", "", opt.merge_dedup, "merge redundancy filter: bloom, sort, or none");
  desc.AddOption("contig_edge_multiplicity", "", opt.contig_edge_multiplicity,
                 "multiplicity assigned to new contig-derived merge edges");
  desc.AddOption("reduce_dedup", "", opt.reduce_dedup, "exact-sort surviving reduced edges, retaining maximum multiplicity");
  desc.AddOption("read_counted_shards", "", opt.read_counted_shards,
                 "prefer records from this many leading read-counted shards during reduction deduplication");
  desc.AddOption("sort_memory_mb", "", opt.sort_memory_mb, "exact-sort record buffer budget in MiB");
  desc.Parse(argc, argv);

  if (opt.long_prefix.empty() || (opt.mode != "index-bench" && opt.output_prefix.empty()) || opt.long_k == 0 ||
      (opt.mode == "reduce" && (opt.short_prefix.empty() || opt.short_k == 0)) ||
      (opt.mode == "index-bench" && opt.short_k == 0)) {
    std::cerr << "Usage: " << argv[0] << " " << std::string(desc) << std::endl;
    std::exit(1);
  }
  if (opt.num_threads <= 0) opt.num_threads = 1;
  if (opt.max_hashes <= 0) opt.max_hashes = 1;
  if (opt.merge_index != "fm" && opt.merge_index != "hash")
    xfatal("Unknown merge index; choose fm or hash\n");
  if ((opt.merge_dedup != "bloom" && opt.merge_dedup != "sort" && opt.merge_dedup != "none") ||
      !std::isfinite(opt.sort_memory_mb) || opt.sort_memory_mb < 0.0625 || opt.sort_memory_mb > 1048576)
    xfatal("Invalid merge dedup backend or sort memory (minimum 0.0625 MiB)\n");
  if (opt.mode == "merge" && opt.merge_dedup == "sort" && opt.human_readable)
    xfatal("Exact sort currently requires --human_readable 0; per-candidate provenance reports are not supported\n");
  if (opt.contig_edge_multiplicity == 0 || opt.contig_edge_multiplicity > kMaxMul)
    xfatal("Invalid contig edge multiplicity; expected 1 to {}\n", kMaxMul);
  if (!std::isfinite(opt.max_bloom_gb) || opt.max_bloom_gb < 0 || opt.max_bloom_gb > 1048576 ||
      !std::isfinite(opt.reserve_gb) || opt.reserve_gb < 0 || opt.reserve_gb > 1048576 ||
      !std::isfinite(opt.host_mem) || opt.host_mem < 0 || opt.host_mem >= 9e18 ||
      !std::isfinite(opt.reduce_fpr) || opt.reduce_fpr <= 0 || opt.reduce_fpr >= 1 ||
      opt.max_hashes > 64 ||
      opt.short_k > kMaxK || opt.long_k > kMaxK)
    xfatal("Invalid Bloom memory, false-positive target, or k-mer size\n");
  opt.human_readable = opt.human_readable != 0 ? 1 : 0;
  return opt;
}

int RunReduce(const Options &opt) {
  using Clock = std::chrono::steady_clock;
  const auto started = Clock::now();
  EdgeIoMetadata short_meta = LoadMetadata(opt.short_prefix);
  EdgeIoMetadata long_meta = LoadMetadata(opt.long_prefix);
  if (short_meta.kmer_size != opt.short_k || long_meta.kmer_size != opt.long_k) {
    xfatal("Invalid k sizes: short edge k {}, long edge k {}, requested short {}, long {}\n",
           short_meta.kmer_size, long_meta.kmer_size, opt.short_k, opt.long_k);
  }

  const unsigned short_edge_len = EdgeLength(short_meta);
  const unsigned long_edge_len = EdgeLength(long_meta);
  if (long_edge_len < short_edge_len) {
    xfatal("Long edge length {} is shorter than short edge length {}\n",
           long_edge_len, short_edge_len);
  }

  auto long_ranges = ReductionRanges(opt.long_prefix, long_meta, opt.num_threads);
  auto short_ranges = ReductionRanges(opt.short_prefix, short_meta, opt.num_threads);
  if (opt.read_counted_shards && (!opt.reduce_dedup || opt.read_counted_shards >= short_meta.num_files))
    throw std::runtime_error("Invalid read-counted shard boundary for reduction deduplication");

  uint64_t long_records = TotalRecords(long_meta);
  uint64_t bloom_insertions = long_records * 2;
  uint64_t bloom_byte_cap = SelectBloomByteCap(opt, 2 * kEdgeBufferBytes * opt.num_threads);
  ReductionBloomShape shape = SelectReductionBloom(bloom_insertions, bloom_byte_cap, opt);
  BloomFilter longer_prefixes(shape.bytes * 8, shape.hashes,
                              (short_edge_len + 31) / 32, opt.num_threads);
  std::ostringstream fpr, requested_fpr;
  fpr << std::scientific << std::setprecision(3) << shape.estimated_fpr;
  requested_fpr << std::scientific << std::setprecision(3) << opt.reduce_fpr;
  xinfo("Reduction Bloom: {} bytes, {} probes, {} active words, estimated FPR {s}, cap {} bytes\n",
        shape.bytes, shape.hashes, (short_edge_len + 31) / 32,
        fpr.str().c_str(), bloom_byte_cap);
  xinfo("Reduction work: {} insertion chunks, {} filtering chunks, {} workers\n",
        long_ranges.size(), short_ranges.size(), opt.num_threads);
  if (shape.estimated_fpr > opt.reduce_fpr)
    xinfo("WARNING: Bloom memory/hash cap prevents requested reduction FPR {s}; increase memory or --max_hashes\n",
          requested_fpr.str().c_str());
  const auto initialized = Clock::now();
  using OriginMap = std::unordered_map<KmerKey, std::vector<LongEdgeOrigin>, KmerKeyHash>;
  // Provenance is diagnostic only; reduction itself uses the packed Bloom keys.
  std::unique_ptr<OriginMap> longer_origin;
  if (opt.human_readable) longer_origin.reset(new OriginMap);

#pragma omp parallel for schedule(dynamic, 1) num_threads(opt.num_threads)
  for (size_t i = 0; i < long_ranges.size(); ++i) {
    const auto &range = long_ranges[i];
    ForEachEdgeRange(opt.long_prefix, long_meta, range.file_id, range.first, range.count,
                     [&](const uint32_t *edge) {
      KmerKey prefix_key = RawKeyFromEdge(edge, 0, short_edge_len);
      KmerKey rc_prefix_key = RawReverseComplementKeyFromEdge(
          edge, long_edge_len - short_edge_len, short_edge_len);
      longer_prefixes.Add(prefix_key);
      longer_prefixes.Add(rc_prefix_key);
      if (longer_origin) {
        std::string long_edge = BasesFromEdge(edge, 0, long_edge_len);
        std::string rc_long_edge = ReverseComplement(long_edge);
#pragma omp critical(longer_origin_insert)
        {
          (*longer_origin)[prefix_key].push_back(LongEdgeOrigin{long_edge, "prefix"});
          (*longer_origin)[rc_prefix_key].push_back(
              LongEdgeOrigin{rc_long_edge, "reverse_complement_prefix"});
        }
      }
    });
  }

  const auto populated = Clock::now();

  std::vector<std::ostringstream> removed_buffers(opt.human_readable ? opt.num_threads : 0);
  std::vector<std::ostringstream> non_removed_buffers(opt.human_readable ? opt.num_threads : 0);
  EdgeOutputWriter writer(opt.output_prefix, opt.short_k, short_meta.words_per_edge,
                          opt.num_threads, opt.human_readable != 0);
  std::unique_ptr<ExactMergeDedup> survivors;
  std::vector<uint64_t> surviving_records(opt.num_threads, 0);
  if (opt.reduce_dedup) {
    if (opt.human_readable) throw std::runtime_error("Reduced-edge dedup requires --human_readable 0");
    uint64_t budget = static_cast<uint64_t>(opt.sort_memory_mb * (1ULL << 20));
    uint64_t available = AvailableMemoryBytes();
    if (opt.host_mem) available = available ? std::min(available, static_cast<uint64_t>(opt.host_mem))
                                          : static_cast<uint64_t>(opt.host_mem);
    // The reduction Bloom remains resident while surviving records are sorted.
    if (available) budget = std::min(budget, available > shape.bytes ? (available - shape.bytes) / 4 : 0);
    survivors = MakeExactMergeDedup((short_edge_len + 31) / 32, opt.num_threads, budget, opt.output_prefix);
    xinfo("Reduced-edge exact dedup: {} workers, {} buffer bytes; {s}\n", opt.num_threads, budget,
          opt.read_counted_shards ? "prefer read-counted records" : "maximum multiplicity");
  }
  constexpr uint64_t source_stride = uint64_t{kMaxMul} + 1;
  auto keep_edge = [&](const uint32_t *edge, const KmerKey &forward, int tid, unsigned file_id) {
    if (!survivors) { writer.Write(edge, tid, "non_removed"); return; }
    ++surviving_records[tid];
    KmerKey reverse = PackedReverseComplement(forward, short_edge_len);
    const auto &key = PackedKeyLess(reverse, forward) ? reverse : forward;
    // The sorter keeps the lowest priority value for each key.
    uint64_t priority = kMaxMul - Multiplicity(edge, short_meta.words_per_edge);
    if (opt.read_counted_shards && file_id >= opt.read_counted_shards) priority += source_stride;
    survivors->AddRecord(key.words, priority, tid);
  };
  std::atomic<uint64_t> kept{0};
  std::atomic<uint64_t> removed{0};
  std::atomic<uint64_t> forward_only{0};
  std::atomic<uint64_t> reverse_only{0};
  std::atomic<uint64_t> no_orientation_hit{0};

#pragma omp parallel for schedule(dynamic, 1) num_threads(opt.num_threads)
  for (size_t i = 0; i < short_ranges.size(); ++i) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    const auto &range = short_ranges[i];
    ForEachEdgeRange(opt.short_prefix, short_meta, range.file_id, range.first, range.count,
                     [&](const uint32_t *edge) {
      KmerKey forward_key = RawKeyFromEdge(edge, 0, short_edge_len);
      bool forward_hit = longer_prefixes.Contains(forward_key);
      if (!forward_hit && !opt.human_readable) {
        keep_edge(edge, forward_key, tid, range.file_id);
        return;
      }
      KmerKey reverse_key = PackedReverseComplement(forward_key, short_edge_len);
      bool reverse_hit = longer_prefixes.Contains(reverse_key);
      if (forward_hit && reverse_hit) {
        if (!opt.human_readable) return;
        ++removed;
        auto forward_origin_it = longer_origin->find(forward_key);
        auto reverse_origin_it = longer_origin->find(reverse_key);
        bool wrote_origin = false;
        if (forward_origin_it != longer_origin->end() &&
            !forward_origin_it->second.empty()) {
          for (const auto &origin : forward_origin_it->second) {
            if (opt.human_readable) removed_buffers[tid] << BasesFromEdge(edge, 0, short_edge_len) << '\t'
                                 << Multiplicity(edge, short_meta.words_per_edge) << '\t'
                                 << "forward_represented_by_longer_edge_"
                                 << origin.side << '\t'
                                 << origin.edge << '\n';
            wrote_origin = true;
          }
        }
        if (reverse_origin_it != longer_origin->end() &&
            !reverse_origin_it->second.empty()) {
          for (const auto &origin : reverse_origin_it->second) {
            if (opt.human_readable) removed_buffers[tid] << BasesFromEdge(edge, 0, short_edge_len) << '\t'
                                 << Multiplicity(edge, short_meta.words_per_edge) << '\t'
                                 << "reverse_complement_represented_by_longer_edge_"
                                 << origin.side << '\t'
                                 << origin.edge << '\n';
            wrote_origin = true;
          }
        }
        if (!wrote_origin) {
          if (opt.human_readable) removed_buffers[tid] << BasesFromEdge(edge, 0, short_edge_len) << '\t'
                               << Multiplicity(edge, short_meta.words_per_edge) << '\t'
                               << "both_orientations_represented_by_longer_edge" << '\t'
                               << "NA" << '\n';
        }
      } else {
        keep_edge(edge, forward_key, tid, range.file_id);
        if (!opt.human_readable) return;
        ++kept;
        const char *keep_reason = "no_orientation_hit";
        if (forward_hit) {
          ++forward_only;
          keep_reason = "forward_only_hit";
        } else if (reverse_hit) {
          ++reverse_only;
          keep_reason = "reverse_complement_only_hit";
        } else {
          ++no_orientation_hit;
        }
        if (opt.human_readable) non_removed_buffers[tid] << BasesFromEdge(edge, 0, short_edge_len) << '\t'
                                 << Multiplicity(edge, short_meta.words_per_edge) << '\t'
                                 << keep_reason << '\n';
      }
    });
  }
  if (survivors) {
    const auto sort_started = Clock::now();
    std::vector<uint64_t> unique(opt.num_threads, 0);
    std::vector<std::vector<uint32_t>> packed_edges(opt.num_threads);
    survivors->FinishUnordered([&](const uint64_t *words, uint64_t priority, int worker) {
      KmerKey key;
      std::copy_n(words, (short_edge_len + 31) / 32, key.words);
      auto &packed = packed_edges[worker];
      PackKeyAsEdge(key, short_edge_len, short_meta.words_per_edge, kMaxMul - priority % source_stride, &packed);
      writer.Write(packed.data(), worker, "non_removed");
      ++unique[worker];
    });
    uint64_t total = 0;
    for (auto count : unique) total += count;
    uint64_t input = 0;
    for (auto count : surviving_records) input += count;
    auto stats = survivors->Stats();
    xinfo("Reduced-edge exact dedup: {} surviving records, {} unique surviving edges, {} duplicates removed, {} temporary bytes written, finalize seconds {}\n",
          input, total, input - total, stats.bytes_written,
          std::chrono::duration<double>(Clock::now() - sort_started).count());
  }
  writer.Finalize();
  const auto filtered = Clock::now();
  xinfo("Reduction seconds: initialize {}, insert {}, filter/write {}, total {}\n",
        std::chrono::duration<double>(initialized - started).count(),
        std::chrono::duration<double>(populated - initialized).count(),
        std::chrono::duration<double>(filtered - populated).count(),
        std::chrono::duration<double>(filtered - started).count());

  if (opt.human_readable) {
    std::ofstream removed_out(opt.output_prefix + ".removed.txt");
    removed_out << "#edge\tmultiplicity\treason\tmatched_long_edge\n";
    for (auto &buf : removed_buffers) removed_out << buf.str();

    std::ofstream non_removed_out(opt.output_prefix + ".non_removed.txt");
    non_removed_out << "#edge\tmultiplicity\treason\n";
    for (auto &buf : non_removed_buffers) non_removed_out << buf.str();

    std::ofstream summary(opt.output_prefix + ".reduction.summary.txt");
    summary << "mode\treduce\n";
  summary << "short_k\t" << opt.short_k << '\n';
  summary << "long_k\t" << opt.long_k << '\n';
  summary << "short_edge_prefix\t" << opt.short_prefix << '\n';
  summary << "long_edge_prefix\t" << opt.long_prefix << '\n';
  summary << "output_edge_prefix\t" << opt.output_prefix << '\n';
  summary << "short_edge_records\t" << TotalRecords(short_meta) << '\n';
  summary << "long_edge_records\t" << long_records << '\n';
  summary << "removed_short_edges\t" << removed.load() << '\n';
  summary << "non_removed_short_edges\t" << kept.load() << '\n';
  summary << "non_removed_forward_only_hits\t" << forward_only.load() << '\n';
  summary << "non_removed_reverse_complement_only_hits\t" << reverse_only.load() << '\n';
  summary << "non_removed_no_orientation_hits\t" << no_orientation_hit.load() << '\n';
  summary << "bloom_insertions\t" << bloom_insertions << '\n';
  summary << "bloom_bits\t" << longer_prefixes.num_bits() << '\n';
  summary << "bloom_bytes\t" << longer_prefixes.num_bytes() << '\n';
  summary << "bloom_hashes\t" << longer_prefixes.num_hashes() << '\n';
  summary << "estimated_false_positive_rate\t"
          << std::setprecision(12)
          << EstimateFpr(longer_prefixes.num_bits(), bloom_insertions,
                         longer_prefixes.num_hashes()) << '\n';
    summary << "representation\tbinary_short_edges_reduced_by_raw_long_edge_and_reverse_complement_prefixes_both_short_orientations_required\n";
  }
  return 0;
}

int RunIndexBenchmark(const Options &opt) {
  using Clock = std::chrono::steady_clock;
  auto started = Clock::now();
  EdgeIoMetadata meta = LoadMetadata(opt.long_prefix);
  unsigned edge_len = EdgeLength(meta);
  if (meta.kmer_size != opt.long_k || opt.short_k >= edge_len)
    xfatal("Invalid k sizes for index benchmark\n");
  MergeQueryIndex prefixes(opt.short_k, opt.merge_index), suffixes(opt.short_k, opt.merge_index);
  std::ifstream contigs(opt.contig_file);
  if (!contigs) xfatal("Cannot read contig-end benchmark input\n");
  std::string name, seq;
  while (ReadNextFasta(&contigs, &name, &seq)) {
    if (seq.size() < opt.short_k) continue;
    PackedMergeContig contig(seq, name, false);
    KmerKey prefix, suffix;
    if (!contig.EndKey(0, opt.short_k, &prefix) ||
        !contig.EndKey(seq.size() - opt.short_k, opt.short_k, &suffix)) continue;
    suffixes.Register(prefix); prefixes.Register(suffix);
  }
  prefixes.Build(opt.num_threads); suffixes.Build(opt.num_threads);
  auto indexed = Clock::now();
  struct Counts { uint64_t records = 0, orientations = 0, prefix = 0, suffix = 0; };
  auto ranges = ReductionRanges(opt.long_prefix, meta, opt.num_threads);
  std::vector<Counts> counts(ranges.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(opt.num_threads)
  for (size_t i = 0; i < ranges.size(); ++i) {
    const auto &range = ranges[i];
    auto &count = counts[i];
    auto query = [&](const KmerKey &key) {
      ++count.orientations;
      count.prefix += prefixes.FindId(PackedKeySlice(key, 0, opt.short_k)) != EndFmIndex::Missing();
      count.suffix += suffixes.FindId(PackedKeySlice(key, edge_len - opt.short_k, opt.short_k)) != EndFmIndex::Missing();
    };
    ForEachEdgeRange(opt.long_prefix, meta, range.file_id, range.first, range.count,
                     [&](const uint32_t *edge) {
      ++count.records;
      auto forward = RawKeyFromEdge(edge, 0, edge_len);
      auto reverse = PackedReverseComplement(forward, edge_len);
      query(forward);
      if (!(reverse == forward)) query(reverse);
    });
  }
  auto finished = Clock::now();
  Counts total;
  for (const auto &count : counts) {
    total.records += count.records; total.orientations += count.orientations;
    total.prefix += count.prefix; total.suffix += count.suffix;
  }
  std::cout << std::setprecision(9) << "{\"backend\":\"" << opt.merge_index
            << "\",\"threads\":" << opt.num_threads << ",\"build_seconds\":"
            << std::chrono::duration<double>(indexed - started).count()
            << ",\"lookup_seconds\":" << std::chrono::duration<double>(finished - indexed).count()
            << ",\"records\":" << total.records << ",\"orientations\":" << total.orientations
            << ",\"prefix_hits\":" << total.prefix << ",\"suffix_hits\":" << total.suffix
            << ",\"prefix_keys\":" << prefixes.size() << ",\"suffix_keys\":" << suffixes.size()
            << ",\"fm_bytes\":" << prefixes.fm_bytes() + suffixes.fm_bytes() << "}\n";
  return 0;
}

int RunMerge(const Options &opt) {
  using Clock = std::chrono::steady_clock;
  const auto started = Clock::now();
  EdgeIoMetadata base_meta = LoadMetadata(opt.long_prefix);
  if (base_meta.kmer_size != opt.long_k) {
    xfatal("Invalid merge k size: edge k {}, requested {}\n",
           base_meta.kmer_size, opt.long_k);
  }
  if (opt.short_k == 0) {
    xfatal("Merge mode requires --short_k for overlap extension\n");
  }

  const unsigned edge_len = EdgeLength(base_meta);
  const unsigned overlap_len = opt.short_k;
  const unsigned active_words = (edge_len + 31) / 32;
  if (overlap_len >= edge_len) {
    xfatal("Invalid overlap length {} for long edge length {}\n", overlap_len, edge_len);
  }

  // Keep source provenance in the shard layout without changing edge records.
  const int candidate_shard_start = opt.merge_dedup == "none" ? opt.num_threads : 0;
  EdgeOutputWriter writer(opt.output_prefix, opt.long_k, base_meta.words_per_edge,
                          opt.num_threads + candidate_shard_start, opt.human_readable != 0);
  uint64_t expected_merge_items = std::max<uint64_t>(1, TotalRecords(base_meta));
  std::unique_ptr<BloomFilter> read_count_edges, emitted_extensions;
  std::unique_ptr<ExactMergeDedup> exact;
  if (opt.merge_dedup == "sort") {
    uint64_t budget = static_cast<uint64_t>(opt.sort_memory_mb * (1ULL << 20));
    uint64_t available = AvailableMemoryBytes();
    if (opt.host_mem) available = available ? std::min(available, static_cast<uint64_t>(opt.host_mem))
                                          : static_cast<uint64_t>(opt.host_mem);
    if (available) budget = std::min(budget, available / 4);
    if (budget < (64ULL << 10)) throw std::runtime_error("Insufficient memory for exact-sort buffers");
    exact = MakeExactMergeDedup(active_words, opt.num_threads, budget, opt.output_prefix);
    xinfo("Merge dedup: exact sort, {} record-buffer bytes, no merge Bloom filters\n", budget);
  } else if (opt.merge_dedup == "bloom") {
    uint64_t bloom_byte_cap = SelectBloomByteCap(opt);
    uint64_t merge_bloom_bytes = std::max<uint64_t>(2ULL << 20, bloom_byte_cap / 2);
    uint64_t read_count_target_bytes = std::max<uint64_t>(1ULL << 20, expected_merge_items * 64ULL);
    uint64_t extension_target_bytes = std::max<uint64_t>(1ULL << 20, expected_merge_items * 64ULL);
    uint64_t read_count_bloom_bytes = std::min(merge_bloom_bytes, read_count_target_bytes);
    uint64_t extension_bloom_bytes = std::min(merge_bloom_bytes, extension_target_bytes);
    read_count_edges.reset(new BloomFilter(read_count_bloom_bytes * 8,
                               ChooseHashCount(read_count_bloom_bytes * 8,
                                               expected_merge_items,
                                               opt.max_hashes)));
    emitted_extensions.reset(new BloomFilter(extension_bloom_bytes * 8,
                                 ChooseHashCount(extension_bloom_bytes * 8,
                                                 expected_merge_items,
                                                 opt.max_hashes)));
    xinfo("Merge dedup: Bloom, {} filter bytes\n", read_count_bloom_bytes + extension_bloom_bytes);
  } else {
    xinfo("Merge dedup: deferred until reduction; no merge Bloom filters or sort\n");
  }
  const auto index_started = Clock::now();
  MergeQueryIndex prefix_queries(overlap_len, opt.merge_index), suffix_queries(overlap_len, opt.merge_index);
  std::vector<std::unique_ptr<PackedMergeContig>> short_contigs;
  uint64_t contig_seen = 0;
  if (!opt.contig_file.empty()) {
    std::ifstream contigs(opt.contig_file);
    std::string name, seq;
    while (ReadNextFasta(&contigs, &name, &seq)) {
      ++contig_seen;
      if (seq.size() < overlap_len) continue;
      std::unique_ptr<PackedMergeContig> contig(new PackedMergeContig(seq, name, opt.human_readable));
      KmerKey prefix, suffix;
      if (!contig->EndKey(0, overlap_len, &prefix) ||
          !contig->EndKey(seq.size() - overlap_len, overlap_len, &suffix)) continue;
      suffix_queries.Register(prefix);
      prefix_queries.Register(suffix);
      short_contigs.push_back(std::move(contig));
    }
  }
  const bool overlap_queries = !short_contigs.empty();
  if (overlap_queries) {
    prefix_queries.Build(opt.num_threads);
    suffix_queries.Build(opt.num_threads);
  }
  const auto indexed = Clock::now();
  if (!overlap_queries) {
    xinfo("Merge index: disabled (long-contig-only merge; no short-contig queries)\n");
  } else {
    xinfo("Merge index: {s}, {} prefix keys, {} suffix keys, {} FM bytes, build seconds {}\n",
          opt.merge_index.c_str(), prefix_queries.size(), suffix_queries.size(),
          prefix_queries.fm_bytes() + suffix_queries.fm_bytes(),
          std::chrono::duration<double>(indexed - index_started).count());
  }

  uint64_t base_seen = 0;
  uint64_t base_added = 0;
  uint64_t base_duplicate = 0;

  uint64_t orientations_scanned = 0;
  uint64_t matched_orientations = 0;
  struct FileMatches {
    struct Hit { size_t prefix, suffix; };
    std::vector<uint64_t> words;
    std::vector<Hit> hits;
    uint64_t records{0};
    uint64_t duplicates{0};
    uint64_t orientations{0};
  };
  const MergeQueryIndex &prefix_lookup = prefix_queries;
  const MergeQueryIndex &suffix_lookup = suffix_queries;
  // Diagnostic duplicate counts observe insertion order, so retain a serial scan there.
  int scan_workers = opt.human_readable ? 1 : opt.num_threads;
  bool split_files = !opt.human_readable && (exact || read_count_edges || overlap_queries) &&
                     base_meta.num_files < static_cast<unsigned>(scan_workers);
  std::vector<EdgeRange> scan_ranges;
  if (split_files) scan_ranges = ReductionRanges(opt.long_prefix, base_meta, scan_workers);
  std::vector<FileMatches> file_matches(split_files ? scan_ranges.size() : base_meta.num_files);

  // With fewer input files than workers, copy their base records once in order,
  // then parallelize the CPU-heavy Bloom/search pass over record ranges. This
  // extra input pass preserves the historical output shard and record order.
  if (split_files) {
#pragma omp parallel for schedule(static) num_threads(scan_workers)
    for (unsigned file = 0; file < base_meta.num_files; ++file)
      ForEachEdgeInFile(opt.long_prefix, base_meta, file, [&](const uint32_t *edge) {
        writer.Write(edge, file % opt.num_threads, "read_count_base");
      });
  }

  auto scan = [&](const EdgeRange &range, FileMatches &matches, bool write_base) {
      auto match_orientation = [&](const KmerKey &edge) {
        ++matches.orientations;
        size_t prefix = prefix_lookup.FindId(PackedKeySlice(edge, 0, overlap_len));
        size_t suffix = suffix_lookup.FindId(PackedKeySlice(edge, edge_len - overlap_len, overlap_len));
        if (prefix != EndFmIndex::Missing() || suffix != EndFmIndex::Missing()) {
          matches.words.insert(matches.words.end(), edge.words, edge.words + active_words);
          matches.hits.push_back({prefix, suffix});
        }
      };
      ForEachEdgeRange(opt.long_prefix, base_meta, range.file_id, range.first, range.count,
                       [&](const uint32_t *edge) {
        ++matches.records;
        if (!exact && !read_count_edges && !overlap_queries) {
          if (write_base) writer.Write(edge, range.file_id % opt.num_threads, "read_count_base");
          return;
        }
        KmerKey forward = RawKeyFromEdge(edge, 0, edge_len);
        KmerKey reverse = PackedReverseComplement(forward, edge_len);
        const KmerKey &key = PackedKeyLess(reverse, forward) ? reverse : forward;
        if (exact) {
          int worker = 0;
#ifdef _OPENMP
          worker = omp_get_thread_num();
#endif
          exact->AddBase(key.words, worker);
        } else if (read_count_edges) {
          if (opt.human_readable && read_count_edges->Contains(key)) ++matches.duplicates;
          read_count_edges->Add(key);
        }
        if (write_base) writer.Write(edge, range.file_id % opt.num_threads, "read_count_base");

        if (overlap_queries) {
          match_orientation(forward);
          if (!(reverse == forward)) match_orientation(reverse);
        }
      });
  };
  if (split_files) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(scan_workers)
    for (size_t i = 0; i < scan_ranges.size(); ++i) scan(scan_ranges[i], file_matches[i], false);
  } else {
#pragma omp parallel for schedule(static) num_threads(scan_workers)
    for (int owner = 0; owner < scan_workers; ++owner) {
      // Files mapping to the same output stream have exactly one worker owner.
      for (unsigned file = owner; file < base_meta.num_files; file += scan_workers)
        scan({file, 0, FileRecordCount(opt.long_prefix, base_meta, file)}, file_matches[file], true);
    }
  }
  const auto scanned = Clock::now();
  if (exact) exact->Check();
  // The barrier above completes the read Bloom before any extension checks it.
  for (FileMatches &matches : file_matches) {
    base_seen += matches.records;
    base_duplicate += matches.duplicates;
    orientations_scanned += matches.orientations;
    matched_orientations += matches.hits.size();
    for (size_t i = 0; i < matches.hits.size(); ++i) {
      KmerKey edge;
      std::copy_n(matches.words.data() + i * active_words, active_words, edge.words);
      auto hit = matches.hits[i];
      if (hit.prefix != EndFmIndex::Missing()) prefix_queries.Hits(hit.prefix).push_back(edge);
      if (hit.suffix != EndFmIndex::Missing()) suffix_queries.Hits(hit.suffix).push_back(edge);
    }
    std::vector<uint64_t>().swap(matches.words);
    std::vector<FileMatches::Hit>().swap(matches.hits);
  }
  base_added = base_seen;
  file_matches.clear();
  const auto consolidated = Clock::now();

  uint64_t contig_orientations = 0;
  uint64_t overlap_matches = 0;
  uint64_t extended_windows_seen = 0;
  uint64_t extended_edges_added = 0;
  uint64_t extended_edges_duplicate = 0;
  uint64_t long_contig_seen = 0;
  uint64_t long_contig_orientations = 0;
  uint64_t long_contig_windows_seen = 0;
  uint64_t long_contig_edges_added = 0;
  uint64_t long_contig_edges_duplicate = 0;
  std::ofstream contig_report;
  if (opt.human_readable) {
    contig_report.open(opt.output_prefix + ".contig_edges.txt");
    contig_report << "#edge\tmultiplicity\taction\tcontig\tdirection\tlong_edge\textension_offset\tduplicate_source\n";
  }

  auto emit_extended_windows = [&](const PackedMergeView &extended,
                                   const std::string &contig_name,
                                   const char *direction,
                                   size_t first_window,
                                   size_t last_window) {
    if (extended.length < edge_len) return;
    std::string long_edge;
    if (opt.human_readable) long_edge = extended.edge ? BasesFromKey(*extended.edge, edge_len) : ".";
    std::vector<uint32_t> packed;
    last_window = std::min(last_window, extended.length - edge_len + 1);
    for (size_t i = first_window; i < last_window; ++i) {
      size_t extension_offset = i - first_window + 1;
      KmerKey key;
      if (!CanonicalPackedWindow(extended, i, edge_len, &key)) continue;
      ++extended_windows_seen;
      if (exact) { exact->AddCandidate(key.words, extended_windows_seen); continue; }
      if (!read_count_edges) {
        PackKeyAsEdge(key, edge_len, base_meta.words_per_edge, opt.contig_edge_multiplicity, &packed);
        writer.Write(packed.data(), candidate_shard_start + extended_edges_added++ % opt.num_threads,
                     "contig_overlap_extension");
        continue;
      }
      const char *duplicate_source = nullptr;
      if (read_count_edges->Contains(key)) {
        duplicate_source = "read_count_base_bloom";
      } else if (emitted_extensions->Contains(key)) {
        duplicate_source = "contig_overlap_extension_bloom";
      }
      if (duplicate_source != nullptr) {
        ++extended_edges_duplicate;
        if (opt.human_readable) contig_report << BasesFromKey(key, edge_len) << '\t'
                      << opt.contig_edge_multiplicity << "\tduplicate\t"
                      << contig_name << '\t' << direction << '\t'
                      << long_edge << '\t' << extension_offset << '\t'
                      << duplicate_source << '\n';
        continue;
      }
      emitted_extensions->Add(key);
      PackKeyAsEdge(key, edge_len, base_meta.words_per_edge,
                    opt.contig_edge_multiplicity, &packed);
      writer.Write(packed.data(), extended_edges_added % opt.num_threads,
                   "contig_overlap_extension");
      ++extended_edges_added;
      if (opt.human_readable) contig_report << BasesFromKey(key, edge_len) << '\t'
                    << opt.contig_edge_multiplicity << "\tadded\t"
                    << contig_name << '\t' << direction << '\t'
                    << long_edge << '\t' << extension_offset << "\t.\n";
    }
  };

  for (const auto &contig_ptr : short_contigs) {
    const PackedMergeContig &contig = *contig_ptr;
    size_t length = contig.bases.size();
    ++contig_orientations;

    KmerKey contig_prefix, contig_suffix;
    contig.EndKey(0, overlap_len, &contig_prefix);
    contig.EndKey(length - overlap_len, overlap_len, &contig_suffix);
    auto suffix_hits = suffix_queries.Find(contig_prefix);
    if (suffix_hits) {
      for (const KmerKey &long_edge : *suffix_hits) {
        ++overlap_matches;
        size_t extra_len = length - overlap_len;
        PackedMergeView extended{contig, &long_edge, edge_len, overlap_len, true, edge_len + extra_len};
        emit_extended_windows(extended, contig.name, "contig_prefix_to_long_suffix", 1, extra_len + 1);
      }
    }

    auto prefix_hits = prefix_queries.Find(contig_suffix);
    if (prefix_hits) {
      for (const KmerKey &long_edge : *prefix_hits) {
        ++overlap_matches;
        size_t extra_len = length - overlap_len;
        PackedMergeView extended{contig, &long_edge, extra_len, 0, false, extra_len + edge_len};
        emit_extended_windows(extended, contig.name, "contig_suffix_to_long_prefix", 0, extra_len);
      }
    }
  }

  const uint64_t long_candidate_start = extended_windows_seen;
  if (!opt.long_contig_file.empty()) {
    std::ifstream long_contigs(opt.long_contig_file);
    if (!long_contigs) throw std::runtime_error("Cannot read long contigs: " + opt.long_contig_file);
    std::string name, seq;
    if (!overlap_queries && !opt.human_readable && !read_count_edges) {
      struct InputContig { std::string name, sequence; uint64_t ordinal; };
      uint64_t next_ordinal = extended_windows_seen + 1;
      std::vector<uint64_t> windows(opt.num_threads, 0);
      std::vector<std::exception_ptr> errors(opt.num_threads);
      xinfo("Long-contig extraction: {} workers, bounded FASTA batches\n", opt.num_threads);
      for (;;) {
        std::vector<InputContig> batch;
        size_t bases = 0;
        while (bases < (8ULL << 20) && ReadNextFasta(&long_contigs, &name, &seq)) {
          ++long_contig_seen;
          if (seq.size() < edge_len) continue;
          ++long_contig_orientations;
          bases += seq.size();
          batch.push_back({name, std::move(seq), next_ordinal});
          next_ordinal += batch.back().sequence.size() - edge_len + 1;
        }
        if (batch.empty()) break;
#pragma omp parallel for schedule(dynamic, 1) num_threads(opt.num_threads)
        for (size_t c = 0; c < batch.size(); ++c) {
          int worker = 0;
#ifdef _OPENMP
          worker = omp_get_thread_num();
#endif
          try {
            const auto &input = batch[c];
            PackedMergeContig contig(input.sequence, input.name, false);
            PackedMergeView view{contig, nullptr, 0, 0, false, contig.bases.size()};
            std::vector<uint32_t> packed;
            for (size_t i = 0; i + edge_len <= input.sequence.size(); ++i) {
              KmerKey key;
              if (!CanonicalPackedWindow(view, i, edge_len, &key)) continue;
              ++windows[worker];
              if (exact) exact->AddRecord(key.words, input.ordinal + i, worker);
              else {
                PackKeyAsEdge(key, edge_len, base_meta.words_per_edge, opt.contig_edge_multiplicity, &packed);
                writer.Write(packed.data(), candidate_shard_start + worker, "long_contig_window");
              }
            }
          } catch (...) { errors[worker] = std::current_exception(); }
        }
        for (const auto &error : errors) if (error) std::rethrow_exception(error);
        if (exact) exact->Check();
      }
      for (auto count : windows) long_contig_windows_seen += count;
      extended_windows_seen += long_contig_windows_seen;
      if (!exact) {
        long_contig_edges_added += long_contig_windows_seen;
        extended_edges_added += long_contig_windows_seen;
      }
    } else {
      while (ReadNextFasta(&long_contigs, &name, &seq)) {
        ++long_contig_seen;
        if (seq.size() < edge_len) continue;
        PackedMergeContig contig(seq, name, opt.human_readable);
        ++long_contig_orientations;
        uint64_t before_seen = extended_windows_seen;
        uint64_t before_added = extended_edges_added;
        uint64_t before_duplicate = extended_edges_duplicate;
        PackedMergeView view{contig, nullptr, 0, 0, false, contig.bases.size()};
        emit_extended_windows(view, contig.name, "long_contig_window", 0, seq.size() - edge_len + 1);
        long_contig_windows_seen += extended_windows_seen - before_seen;
        long_contig_edges_added += extended_edges_added - before_added;
        long_contig_edges_duplicate += extended_edges_duplicate - before_duplicate;
      }
    }
  }
  if (exact) {
    const auto generated = Clock::now();
    std::vector<uint32_t> packed;
    exact->Finish([&](const uint64_t *words, uint64_t ordinal) {
      KmerKey key;
      std::copy_n(words, active_words, key.words);
      PackKeyAsEdge(key, edge_len, base_meta.words_per_edge, opt.contig_edge_multiplicity, &packed);
      writer.Write(packed.data(), extended_edges_added % opt.num_threads, "contig_overlap_extension");
      ++extended_edges_added;
      if (ordinal > long_candidate_start) ++long_contig_edges_added;
    });
    extended_edges_duplicate = extended_windows_seen - extended_edges_added;
    long_contig_edges_duplicate = long_contig_windows_seen - long_contig_edges_added;
    auto stats = exact->Stats();
    xinfo("Exact sort: {} candidates, {} new edges, {} duplicate candidates, {} runs, {} temporary bytes written, {} temporary bytes read, finalize seconds {}\n",
          extended_windows_seen, extended_edges_added, extended_edges_duplicate, stats.runs_written,
          stats.bytes_written, stats.bytes_read,
          std::chrono::duration<double>(Clock::now() - generated).count());
  }
  writer.Finalize();
  const auto finished = Clock::now();
  xinfo("Merge seconds: setup {}, index {}, scan {}, consolidate {}, extend/write {}, total {} ({} scan workers)\n",
        std::chrono::duration<double>(index_started - started).count(),
        std::chrono::duration<double>(indexed - index_started).count(),
        std::chrono::duration<double>(scanned - indexed).count(),
        std::chrono::duration<double>(consolidated - scanned).count(),
        std::chrono::duration<double>(finished - consolidated).count(),
        std::chrono::duration<double>(finished - started).count(), scan_workers);

  if (opt.human_readable) {
    std::ofstream summary(opt.output_prefix + ".merge.summary.txt");
    summary << "mode\tmerge\n";
  summary << "long_k\t" << opt.long_k << '\n';
  summary << "short_k\t" << opt.short_k << '\n';
  summary << "overlap_length\t" << overlap_len << '\n';
  summary << "base_edge_prefix\t" << opt.long_prefix << '\n';
  summary << "output_edge_prefix\t" << opt.output_prefix << '\n';
  summary << "base_edges_seen\t" << base_seen << '\n';
  summary << "base_edges_added\t" << base_added << '\n';
  summary << "base_edges_bloom_hits\t" << base_duplicate << '\n';
  summary << "merge_read_count_bloom_bits\t" << read_count_edges->num_bits() << '\n';
  summary << "merge_read_count_bloom_bytes\t" << read_count_edges->num_bytes() << '\n';
  summary << "merge_read_count_bloom_hashes\t" << read_count_edges->num_hashes() << '\n';
  summary << "merge_read_count_bloom_estimated_fpr\t"
          << std::setprecision(12)
          << EstimateFpr(read_count_edges->num_bits(), expected_merge_items,
                         read_count_edges->num_hashes()) << '\n';
  summary << "merge_extension_bloom_bits\t" << emitted_extensions->num_bits() << '\n';
  summary << "merge_extension_bloom_bytes\t" << emitted_extensions->num_bytes() << '\n';
  summary << "merge_extension_bloom_hashes\t" << emitted_extensions->num_hashes() << '\n';
  summary << "merge_extension_bloom_estimated_fpr\t"
          << std::setprecision(12)
          << EstimateFpr(emitted_extensions->num_bits(), expected_merge_items,
                         emitted_extensions->num_hashes()) << '\n';
  summary << "long_oriented_edges_scanned\t" << orientations_scanned << '\n';
  summary << "long_oriented_edges_matched\t" << matched_orientations << '\n';
  summary << "prefix_query_keys\t" << prefix_queries.size() << '\n';
  summary << "suffix_query_keys\t" << suffix_queries.size() << '\n';
  summary << "scan_workers\t" << scan_workers << '\n';
  summary << "contigs_seen\t" << contig_seen << '\n';
  summary << "contig_orientations_scanned\t" << contig_orientations << '\n';
  summary << "contig_orientation_modes\trev\n";
  summary << "long_contigs_seen\t" << long_contig_seen << '\n';
  summary << "long_contig_orientations_scanned\t" << long_contig_orientations << '\n';
  summary << "long_contig_orientation_modes\trev\n";
  summary << "long_contig_windows_seen\t" << long_contig_windows_seen << '\n';
  summary << "long_contig_edges_added\t" << long_contig_edges_added << '\n';
  summary << "long_contig_edges_duplicate\t" << long_contig_edges_duplicate << '\n';
  summary << "overlap_matches\t" << overlap_matches << '\n';
  summary << "extended_windows_seen\t" << extended_windows_seen << '\n';
  summary << "extended_edges_added\t" << extended_edges_added << '\n';
  summary << "extended_edges_duplicate\t" << extended_edges_duplicate << '\n';
  summary << "contig_edge_multiplicity\t" << opt.contig_edge_multiplicity << '\n';
  summary << "output_edges_total\t" << base_added + extended_edges_added << '\n';
    summary << "representation\tbinary_edges_read_count_plus_exact_k_overlap_extensions_bloom_redundancy_filter\n";
  }
  return 0;
}

}  // namespace

int main_prefix_bloom(int argc, char **argv) {
  try {
    Options opt = ParseOptions(argc, argv);
    if (opt.mode == "merge") return RunMerge(opt);
    if (opt.mode == "index-bench") return RunIndexBenchmark(opt);
    if (opt.mode == "reduce") return RunReduce(opt);
    xfatal("Unknown prefixbf mode {s}\n", opt.mode.c_str());
    return 1;
  } catch (const std::exception &error) {
    std::cerr << "prefixbf: " << error.what() << '\n';
    return 1;
  }
}
