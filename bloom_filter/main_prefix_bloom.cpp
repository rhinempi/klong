#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "definitions.h"
#include "utils/utils.h"
#include "sequence/io/edge/edge_io_meta.h"
#include "utils/options_description.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr uint64_t kDefaultMaxBloomBytes = 60ULL * 1024 * 1024 * 1024;
constexpr uint64_t kSplitMixSeed = 0x9e3779b97f4a7c15ULL;
constexpr uint32_t kContigEdgeMultiplicity = 50;

struct Options {
  std::string mode{"reduce"};
  std::string short_prefix;
  std::string long_prefix;
  std::string output_prefix;
  std::string contig_file;
  unsigned short_k{0};
  unsigned long_k{0};
  int num_threads{1};
  double max_bloom_gb{60.0};
  double reserve_gb{1.0};
  int max_hashes{16};
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

uint64_t HashKey(const KmerKey &key, uint64_t seed) {
  uint64_t h = seed;
  for (int i = 0; i < KmerKey::kWords; ++i) {
    h = SplitMix64(h ^ key.words[i]);
  }
  return h;
}

struct KmerKeyHash {
  size_t operator()(const KmerKey &key) const {
    return static_cast<size_t>(HashKey(key, 0x9e3779b97f4a7c15ULL));
  }
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

bool CanonicalKeyFromString(const std::string &seq, unsigned offset,
                            unsigned len, KmerKey *out) {
  KmerKey fwd, rev;
  for (unsigned i = 0; i < len; ++i) {
    uint8_t b;
    switch (seq[offset + i]) {
      case 'A': case 'a': b = 0; break;
      case 'C': case 'c': b = 1; break;
      case 'G': case 'g': b = 2; break;
      case 'T': case 't': b = 3; break;
      default: return false;
    }
    SetPackedBase(&fwd, i, b);
    SetPackedBase(&rev, len - 1 - i, 3 - b);
  }
  *out = std::memcmp(fwd.words, rev.words, sizeof(fwd.words)) <= 0 ? fwd : rev;
  return true;
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
    if (!is) return 0;
    return static_cast<uint64_t>(is.tellg()) /
           (sizeof(uint32_t) * meta.words_per_edge);
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
  while (is >> key >> value_kb >> unit) {
    if (key == "MemAvailable:") return value_kb * 1024ULL;
  }
  return 0;
}

class BloomFilter {
 public:
  BloomFilter(uint64_t num_bits, int num_hashes)
      : num_bits_(std::max<uint64_t>(64, num_bits)),
        num_words_((num_bits_ + 63) / 64),
        num_hashes_(std::max(1, num_hashes)),
        words_(new std::atomic<uint64_t>[num_words_]) {
    for (uint64_t i = 0; i < num_words_; ++i) words_[i].store(0, std::memory_order_relaxed);
  }

  void Add(const KmerKey &key) {
    uint64_t h1 = HashKey(key, 0x243f6a8885a308d3ULL);
    uint64_t h2 = HashKey(key, 0x13198a2e03707344ULL) | 1ULL;
    for (int i = 0; i < num_hashes_; ++i) {
      uint64_t bit = (h1 + static_cast<uint64_t>(i) * h2) % num_bits_;
      words_[bit >> 6].fetch_or(1ULL << (bit & 63), std::memory_order_relaxed);
    }
  }

  bool Contains(const KmerKey &key) const {
    uint64_t h1 = HashKey(key, 0x243f6a8885a308d3ULL);
    uint64_t h2 = HashKey(key, 0x13198a2e03707344ULL) | 1ULL;
    for (int i = 0; i < num_hashes_; ++i) {
      uint64_t bit = (h1 + static_cast<uint64_t>(i) * h2) % num_bits_;
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
  std::unique_ptr<std::atomic<uint64_t>[]> words_;
};

class EdgeOutputWriter {
 public:
  EdgeOutputWriter(const std::string &prefix, uint32_t k, uint32_t words_per_edge,
                   int num_files)
      : prefix_(prefix), k_(k), words_per_edge_(words_per_edge),
        counts_(std::max(1, num_files), 0) {
    int n = std::max(1, num_files);
    for (int i = 0; i < n; ++i) {
      binary_.emplace_back(new std::ofstream(prefix + ".edges." + std::to_string(i),
                                             std::ofstream::binary | std::ofstream::out));
      text_.emplace_back(new std::ofstream(prefix + ".edges." + std::to_string(i) + ".txt"));
      *text_.back() << "#edge\tmultiplicity\tsource\n";
    }
  }

  void Write(const uint32_t *edge, int file_id, const char *source) {
    file_id %= static_cast<int>(binary_.size());
    binary_[file_id]->write(reinterpret_cast<const char *>(edge),
                            sizeof(uint32_t) * words_per_edge_);
    *text_[file_id] << BasesFromEdge(edge, 0, k_ + 1) << '\t'
                    << Multiplicity(edge, words_per_edge_) << '\t' << source << '\n';
    ++counts_[file_id];
  }

  void Finalize() {
    for (auto &f : binary_) f->close();
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
  }

 private:
  std::string prefix_;
  uint32_t k_;
  uint32_t words_per_edge_;
  std::vector<std::unique_ptr<std::ofstream>> binary_;
  std::vector<std::unique_ptr<std::ofstream>> text_;
  std::vector<int64_t> counts_;
};

void ForEachEdgeInFile(const std::string &prefix, const EdgeIoMetadata &meta,
                       unsigned file_id,
                       const std::function<void(const uint32_t *)> &fn) {
  uint64_t n = FileRecordCount(prefix, meta, file_id);
  if (n == 0) return;

  std::ifstream is(prefix + ".edges." + std::to_string(file_id),
                   std::ifstream::binary | std::ifstream::in);
  if (!is) xfatal("Cannot open {s}.edges.{}\n", prefix.c_str(), file_id);

  std::vector<uint32_t> edge(meta.words_per_edge);
  for (uint64_t i = 0; i < n; ++i) {
    is.read(reinterpret_cast<char *>(edge.data()),
            sizeof(uint32_t) * meta.words_per_edge);
    if (!is) xfatal("Unexpected EOF reading {s}.edges.{}\n", prefix.c_str(), file_id);
    fn(edge.data());
  }
}

int ChooseHashCount(uint64_t bits, uint64_t items, int max_hashes) {
  if (items == 0) return 1;
  double optimal = (static_cast<double>(bits) / static_cast<double>(items)) * std::log(2.0);
  return std::max(1, std::min(max_hashes, static_cast<int>(std::round(optimal))));
}

double EstimateFpr(uint64_t bits, uint64_t items, int hashes) {
  if (bits == 0) return 1.0;
  double exponent = -static_cast<double>(hashes) * static_cast<double>(items) /
                    static_cast<double>(bits);
  return std::pow(1.0 - std::exp(exponent), hashes);
}

uint64_t SelectBloomByteCap(const Options &opt) {
  uint64_t requested = static_cast<uint64_t>(opt.max_bloom_gb * 1024.0 * 1024.0 * 1024.0);
  uint64_t reserve = static_cast<uint64_t>(opt.reserve_gb * 1024.0 * 1024.0 * 1024.0);
  uint64_t available = AvailableMemoryBytes();
  if (requested == 0) requested = kDefaultMaxBloomBytes;
  if (available == 0) return requested;
  if (available <= reserve + (256ULL << 20)) return 256ULL << 20;
  return std::min(requested, available - reserve);
}

bool ReadNextFasta(std::ifstream *in, std::string *name, std::string *seq) {
  seq->clear();
  std::string line;
  if (name->empty()) {
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
      *name = line.substr(1);
      return true;
    }
    seq->append(line);
  }
  name->clear();
  return !seq->empty();
}

Options ParseOptions(int argc, char **argv) {
  Options opt;
  OptionsDescription desc;
  desc.AddOption("mode", "M", opt.mode, "reduce or merge");
  desc.AddOption("short_prefix", "s", opt.short_prefix, "prefix of shorter k binary edge files");
  desc.AddOption("long_prefix", "l", opt.long_prefix, "prefix of longer/current k binary edge files");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix for edge files");
  desc.AddOption("contig", "c", opt.contig_file, "contigs used by merge mode");
  desc.AddOption("short_k", "k", opt.short_k, "shorter k-mer size");
  desc.AddOption("long_k", "K", opt.long_k, "longer/current k-mer size");
  desc.AddOption("num_cpu_threads", "t", opt.num_threads, "number of CPU threads/output files");
  desc.AddOption("max_bloom_gb", "m", opt.max_bloom_gb, "max Bloom filter memory in GB");
  desc.AddOption("reserve_gb", "r", opt.reserve_gb, "RAM to leave unused in GB");
  desc.AddOption("max_hashes", "H", opt.max_hashes, "maximum number of Bloom hash probes");
  desc.Parse(argc, argv);

  if (opt.long_prefix.empty() || opt.output_prefix.empty() || opt.long_k == 0 ||
      (opt.mode == "reduce" && (opt.short_prefix.empty() || opt.short_k == 0))) {
    std::cerr << "Usage: " << argv[0] << " " << std::string(desc) << std::endl;
    std::exit(1);
  }
  if (opt.num_threads <= 0) opt.num_threads = 1;
  if (opt.max_hashes <= 0) opt.max_hashes = 1;
  return opt;
}

int RunReduce(const Options &opt) {
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

  uint64_t long_records = TotalRecords(long_meta);
  uint64_t bloom_insertions = long_records * 2;
  uint64_t bloom_byte_cap = SelectBloomByteCap(opt);
  uint64_t target_bytes = std::max<uint64_t>(64ULL << 20, bloom_insertions * 64ULL);
  uint64_t bloom_bytes = std::min(bloom_byte_cap, target_bytes);
  BloomFilter longer_prefixes(bloom_bytes * 8,
                              ChooseHashCount(bloom_bytes * 8,
                                              std::max<uint64_t>(1, bloom_insertions),
                                              opt.max_hashes));

#pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(long_meta.num_files); ++file_id) {
    ForEachEdgeInFile(opt.long_prefix, long_meta, file_id, [&](const uint32_t *edge) {
      longer_prefixes.Add(CanonicalKeyFromEdge(edge, 0, short_edge_len));
      longer_prefixes.Add(CanonicalKeyFromEdge(edge, long_edge_len - short_edge_len,
                                               short_edge_len));
    });
  }

  std::vector<std::ostringstream> removed_buffers(opt.num_threads);
  std::vector<std::ostringstream> non_removed_buffers(opt.num_threads);
  EdgeOutputWriter writer(opt.output_prefix, opt.short_k, short_meta.words_per_edge,
                          opt.num_threads);
  std::atomic<uint64_t> kept{0};
  std::atomic<uint64_t> removed{0};

#pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(short_meta.num_files); ++file_id) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::vector<std::vector<uint32_t>> local_kept;
    ForEachEdgeInFile(opt.short_prefix, short_meta, file_id, [&](const uint32_t *edge) {
      KmerKey key = CanonicalKeyFromEdge(edge, 0, short_edge_len);
      if (longer_prefixes.Contains(key)) {
        ++removed;
        removed_buffers[tid] << BasesFromEdge(edge, 0, short_edge_len) << '\t'
                             << Multiplicity(edge, short_meta.words_per_edge) << '\t'
                             << "represented_by_longer_edge_prefix" << '\n';
      } else {
        ++kept;
        local_kept.push_back(std::vector<uint32_t>(edge, edge + short_meta.words_per_edge));
        non_removed_buffers[tid] << BasesFromEdge(edge, 0, short_edge_len) << '\t'
                                 << Multiplicity(edge, short_meta.words_per_edge) << '\t'
                                 << "non_removed" << '\n';
      }
    });
    for (const auto &edge : local_kept) writer.Write(edge.data(), tid, "non_removed");
  }
  writer.Finalize();

  std::ofstream removed_out(opt.output_prefix + ".removed.txt");
  removed_out << "#edge\tmultiplicity\treason\n";
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
  summary << "bloom_insertions\t" << bloom_insertions << '\n';
  summary << "bloom_bits\t" << longer_prefixes.num_bits() << '\n';
  summary << "bloom_bytes\t" << longer_prefixes.num_bytes() << '\n';
  summary << "bloom_hashes\t" << longer_prefixes.num_hashes() << '\n';
  summary << "estimated_false_positive_rate\t"
          << std::setprecision(12)
          << EstimateFpr(longer_prefixes.num_bits(), bloom_insertions,
                         longer_prefixes.num_hashes()) << '\n';
  summary << "representation\tbinary_short_edges_reduced_by_long_edge_prefixes\n";
  return 0;
}

int RunMerge(const Options &opt) {
  EdgeIoMetadata base_meta = LoadMetadata(opt.long_prefix);
  if (base_meta.kmer_size != opt.long_k) {
    xfatal("Invalid merge k size: edge k {}, requested {}\n",
           base_meta.kmer_size, opt.long_k);
  }
  const unsigned edge_len = EdgeLength(base_meta);
  EdgeOutputWriter writer(opt.output_prefix, opt.long_k, base_meta.words_per_edge,
                          opt.num_threads);
  std::unordered_set<KmerKey, KmerKeyHash> emitted;
  uint64_t base_seen = 0;
  uint64_t base_added = 0;
  uint64_t base_duplicate = 0;

  for (int file_id = 0; file_id < static_cast<int>(base_meta.num_files); ++file_id) {
    ForEachEdgeInFile(opt.long_prefix, base_meta, file_id, [&](const uint32_t *edge) {
      ++base_seen;
      KmerKey key = CanonicalKeyFromEdge(edge, 0, edge_len);
      if (!emitted.insert(key).second) {
        ++base_duplicate;
        return;
      }
      writer.Write(edge, file_id, "read_count_base");
      ++base_added;
    });
  }

  uint64_t contig_seen = 0;
  uint64_t contig_added = 0;
  uint64_t contig_duplicate = 0;
  std::ofstream contig_report(opt.output_prefix + ".contig_edges.txt");
  contig_report << "#edge\tmultiplicity\taction\n";
  if (!opt.contig_file.empty()) {
    std::ifstream contigs(opt.contig_file);
    std::string name, seq;
    std::vector<uint32_t> packed;
    while (ReadNextFasta(&contigs, &name, &seq)) {
      if (seq.size() < edge_len) continue;
      for (size_t i = 0; i + edge_len <= seq.size(); ++i) {
        KmerKey key;
        if (!CanonicalKeyFromString(seq, i, edge_len, &key)) continue;
        ++contig_seen;
        if (!emitted.insert(key).second) {
          ++contig_duplicate;
          contig_report << BasesFromKey(key, edge_len) << '\t'
                        << kContigEdgeMultiplicity << "\tduplicate\n";
          continue;
        }
        PackKeyAsEdge(key, edge_len, base_meta.words_per_edge,
                      kContigEdgeMultiplicity, &packed);
        writer.Write(packed.data(), contig_added % opt.num_threads,
                     "previous_contig_window");
        ++contig_added;
        contig_report << BasesFromKey(key, edge_len) << '\t'
                      << kContigEdgeMultiplicity << "\tadded\n";
      }
    }
  }
  writer.Finalize();

  std::ofstream summary(opt.output_prefix + ".merge.summary.txt");
  summary << "mode\tmerge\n";
  summary << "long_k\t" << opt.long_k << '\n';
  summary << "base_edge_prefix\t" << opt.long_prefix << '\n';
  summary << "output_edge_prefix\t" << opt.output_prefix << '\n';
  summary << "base_edges_seen\t" << base_seen << '\n';
  summary << "base_edges_added\t" << base_added << '\n';
  summary << "base_edges_duplicate\t" << base_duplicate << '\n';
  summary << "contig_edges_seen\t" << contig_seen << '\n';
  summary << "contig_edges_added\t" << contig_added << '\n';
  summary << "contig_edges_duplicate\t" << contig_duplicate << '\n';
  summary << "contig_edge_multiplicity\t" << kContigEdgeMultiplicity << '\n';
  summary << "output_edges_total\t" << emitted.size() << '\n';
  summary << "representation\tbinary_edges_read_count_plus_previous_contig_windows\n";
  return 0;
}

}  // namespace

int main_prefix_bloom(int argc, char **argv) {
  Options opt = ParseOptions(argc, argv);
  if (opt.mode == "merge") return RunMerge(opt);
  if (opt.mode == "reduce") return RunReduce(opt);
  xfatal("Unknown prefixbf mode {s}\n", opt.mode.c_str());
  return 1;
}
