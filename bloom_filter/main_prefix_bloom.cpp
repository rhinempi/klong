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
#include <unordered_map>
#include <vector>

#include "definitions.h"
#include "kmlib/kmcompactvector.h"
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
  std::string long_contig_file;
  unsigned short_k{0};
  unsigned long_k{0};
  int num_threads{1};
  double max_bloom_gb{60.0};
  double reserve_gb{1.0};
  int max_hashes{16};
  int human_readable{0};
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

KmerKey PackedReverseComplement(const KmerKey &key, unsigned length) {
  KmerKey result;
  if (length <= 32) {
    uint64_t x = ~key.words[0];
    x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
    x = ((x >> 4) & 0x0f0f0f0f0f0f0f0fULL) | ((x & 0x0f0f0f0f0f0f0f0fULL) << 4);
    x = ((x >> 8) & 0x00ff00ff00ff00ffULL) | ((x & 0x00ff00ff00ff00ffULL) << 8);
    x = ((x >> 16) & 0x0000ffff0000ffffULL) | ((x & 0x0000ffff0000ffffULL) << 16);
    x = (x >> 32) | (x << 32);
    result.words[0] = length == 32 ? x : x << (64 - 2 * length);
  } else {
    for (unsigned i = 0; i < length; ++i) {
      SetPackedBase(&result, i, 3 - PackedKeyBase(key, length - 1 - i));
    }
  }
  return result;
}

class MergeQueryIndex {
 public:
  explicit MergeQueryIndex(unsigned length) : compact_(length <= 32) {}

  void Register(const KmerKey &key) {
    if (compact_) small_.emplace(key.words[0], std::vector<KmerKey>{});
    else large_.emplace(key, std::vector<KmerKey>{});
  }

  std::vector<KmerKey> *Find(const KmerKey &key) {
    if (compact_) {
      auto it = small_.find(key.words[0]);
      return it == small_.end() ? nullptr : &it->second;
    }
    auto it = large_.find(key);
    return it == large_.end() ? nullptr : &it->second;
  }

  bool Has(const KmerKey &key) const {
    return compact_ ? small_.find(key.words[0]) != small_.end()
                    : large_.find(key) != large_.end();
  }

  size_t size() const { return compact_ ? small_.size() : large_.size(); }

 private:
  bool compact_;
  std::unordered_map<uint64_t, std::vector<KmerKey>> small_;
  std::unordered_map<KmerKey, std::vector<KmerKey>, KmerKeyHash> large_;
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
  for (unsigned i = 0; i < len; ++i) {
    SetPackedBase(&key, i, GetBase(edge, offset + i));
  }
  return key;
}

KmerKey RawReverseComplementKeyFromEdge(const uint32_t *edge, unsigned offset,
                                        unsigned len) {
  KmerKey key;
  for (unsigned i = 0; i < len; ++i) {
    SetPackedBase(&key, i,
                  static_cast<uint8_t>(3 - GetBase(edge, offset + len - 1 - i)));
  }
  return key;
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
                   int num_files, bool human_readable)
      : prefix_(prefix), k_(k), words_per_edge_(words_per_edge),
        human_readable_(human_readable), counts_(std::max(1, num_files), 0) {
    int n = std::max(1, num_files);
    for (int i = 0; i < n; ++i) {
      binary_.emplace_back(new std::ofstream(prefix + ".edges." + std::to_string(i),
                                             std::ofstream::binary | std::ofstream::out));
      if (human_readable_) {
        text_.emplace_back(new std::ofstream(prefix + ".edges." + std::to_string(i) + ".txt"));
        *text_.back() << "#edge\tmultiplicity\tsource\n";
      }
    }
  }

  void Write(const uint32_t *edge, int file_id, const char *source) {
    file_id %= static_cast<int>(binary_.size());
    binary_[file_id]->write(reinterpret_cast<const char *>(edge),
                            sizeof(uint32_t) * words_per_edge_);
    if (human_readable_) {
      *text_[file_id] << BasesFromEdge(edge, 0, k_ + 1) << '\t'
                      << Multiplicity(edge, words_per_edge_) << '\t' << source << '\n';
    }
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
  bool human_readable_;
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
  desc.AddOption("mode", "M", opt.mode, "reduce or merge");
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
  desc.AddOption("human_readable", "", opt.human_readable, "write human-readable debug reports");
  desc.Parse(argc, argv);

  if (opt.long_prefix.empty() || opt.output_prefix.empty() || opt.long_k == 0 ||
      (opt.mode == "reduce" && (opt.short_prefix.empty() || opt.short_k == 0))) {
    std::cerr << "Usage: " << argv[0] << " " << std::string(desc) << std::endl;
    std::exit(1);
  }
  if (opt.num_threads <= 0) opt.num_threads = 1;
  if (opt.max_hashes <= 0) opt.max_hashes = 1;
  opt.human_readable = opt.human_readable != 0 ? 1 : 0;
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
  using OriginMap = std::unordered_map<KmerKey, std::vector<LongEdgeOrigin>, KmerKeyHash>;
  // Provenance is diagnostic only; reduction itself uses the packed Bloom keys.
  std::unique_ptr<OriginMap> longer_origin;
  if (opt.human_readable) longer_origin.reset(new OriginMap);

#pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(long_meta.num_files); ++file_id) {
    ForEachEdgeInFile(opt.long_prefix, long_meta, file_id, [&](const uint32_t *edge) {
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

  std::vector<std::ostringstream> removed_buffers(opt.human_readable ? opt.num_threads : 0);
  std::vector<std::ostringstream> non_removed_buffers(opt.human_readable ? opt.num_threads : 0);
  EdgeOutputWriter writer(opt.output_prefix, opt.short_k, short_meta.words_per_edge,
                          opt.num_threads, opt.human_readable != 0);
  std::atomic<uint64_t> kept{0};
  std::atomic<uint64_t> removed{0};
  std::atomic<uint64_t> forward_only{0};
  std::atomic<uint64_t> reverse_only{0};
  std::atomic<uint64_t> no_orientation_hit{0};

#pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(short_meta.num_files); ++file_id) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::vector<std::vector<uint32_t>> local_kept;
    ForEachEdgeInFile(opt.short_prefix, short_meta, file_id, [&](const uint32_t *edge) {
      KmerKey forward_key = RawKeyFromEdge(edge, 0, short_edge_len);
      KmerKey reverse_key = RawReverseComplementKeyFromEdge(edge, 0, short_edge_len);
      bool forward_hit = longer_prefixes.Contains(forward_key);
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
        local_kept.push_back(std::vector<uint32_t>(edge, edge + short_meta.words_per_edge));
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
    for (const auto &edge : local_kept) writer.Write(edge.data(), tid, "non_removed");
  }
  writer.Finalize();

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

int RunMerge(const Options &opt) {
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
  if (overlap_len >= edge_len) {
    xfatal("Invalid overlap length {} for long edge length {}\n", overlap_len, edge_len);
  }

  EdgeOutputWriter writer(opt.output_prefix, opt.long_k, base_meta.words_per_edge,
                          opt.num_threads, opt.human_readable != 0);
  uint64_t bloom_byte_cap = SelectBloomByteCap(opt);
  uint64_t merge_bloom_bytes = std::max<uint64_t>(2ULL << 20, bloom_byte_cap / 2);
  uint64_t expected_merge_items = std::max<uint64_t>(1, TotalRecords(base_meta));
  uint64_t read_count_target_bytes = std::max<uint64_t>(1ULL << 20, expected_merge_items * 64ULL);
  uint64_t extension_target_bytes = std::max<uint64_t>(1ULL << 20, expected_merge_items * 64ULL);
  uint64_t read_count_bloom_bytes = std::min(merge_bloom_bytes, read_count_target_bytes);
  uint64_t extension_bloom_bytes = std::min(merge_bloom_bytes, extension_target_bytes);
  BloomFilter read_count_edges(read_count_bloom_bytes * 8,
                               ChooseHashCount(read_count_bloom_bytes * 8,
                                               expected_merge_items,
                                               opt.max_hashes));
  BloomFilter emitted_extensions(extension_bloom_bytes * 8,
                                 ChooseHashCount(extension_bloom_bytes * 8,
                                                 expected_merge_items,
                                                 opt.max_hashes));
  MergeQueryIndex prefix_queries(overlap_len), suffix_queries(overlap_len);
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

  uint64_t base_seen = 0;
  uint64_t base_added = 0;
  uint64_t base_duplicate = 0;

  uint64_t orientations_scanned = 0;
  uint64_t matched_orientations = 0;
  struct FileMatches {
    std::vector<KmerKey> edges;
    uint64_t records{0};
    uint64_t duplicates{0};
    uint64_t orientations{0};
  };
  std::vector<FileMatches> file_matches(base_meta.num_files);
  const MergeQueryIndex &prefix_lookup = prefix_queries;
  const MergeQueryIndex &suffix_lookup = suffix_queries;
  // Diagnostic duplicate counts observe insertion order, so retain a serial scan there.
  int scan_workers = opt.human_readable ? 1 : std::max(1, std::min(
      opt.num_threads, static_cast<int>(base_meta.num_files)));

#pragma omp parallel for schedule(static) num_threads(scan_workers)
  for (int owner = 0; owner < scan_workers; ++owner) {
    // Files mapping to the same output stream have exactly one worker owner.
    for (int file_id = owner; file_id < static_cast<int>(base_meta.num_files);
         file_id += scan_workers) {
      FileMatches &matches = file_matches[file_id];
      auto match_orientation = [&](const KmerKey &edge) {
        ++matches.orientations;
        if (prefix_lookup.Has(PackedKeySlice(edge, 0, overlap_len)) ||
            suffix_lookup.Has(PackedKeySlice(edge, edge_len - overlap_len, overlap_len))) {
          matches.edges.push_back(edge);
        }
      };
      ForEachEdgeInFile(opt.long_prefix, base_meta, file_id, [&](const uint32_t *edge) {
        ++matches.records;
        KmerKey forward = RawKeyFromEdge(edge, 0, edge_len);
        KmerKey reverse = PackedReverseComplement(forward, edge_len);
        const KmerKey &key = PackedKeyLess(reverse, forward) ? reverse : forward;
        if (opt.human_readable && read_count_edges.Contains(key)) ++matches.duplicates;
        read_count_edges.Add(key);
        writer.Write(edge, file_id, "read_count_base");

        match_orientation(forward);
        if (!(reverse == forward)) match_orientation(reverse);
      });
    }
  }
  // The barrier above completes the read Bloom before any extension checks it.
  for (const FileMatches &matches : file_matches) {
    base_seen += matches.records;
    base_duplicate += matches.duplicates;
    orientations_scanned += matches.orientations;
    matched_orientations += matches.edges.size();
    for (const KmerKey &edge : matches.edges) {
      auto prefix_hits = prefix_queries.Find(PackedKeySlice(edge, 0, overlap_len));
      auto suffix_hits = suffix_queries.Find(PackedKeySlice(edge, edge_len - overlap_len, overlap_len));
      if (prefix_hits) prefix_hits->push_back(edge);
      if (suffix_hits) suffix_hits->push_back(edge);
    }
  }
  base_added = base_seen;
  file_matches.clear();

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
      const char *duplicate_source = nullptr;
      if (read_count_edges.Contains(key)) {
        duplicate_source = "read_count_base_bloom";
      } else if (emitted_extensions.Contains(key)) {
        duplicate_source = "contig_overlap_extension_bloom";
      }
      if (duplicate_source != nullptr) {
        ++extended_edges_duplicate;
        if (opt.human_readable) contig_report << BasesFromKey(key, edge_len) << '\t'
                      << kContigEdgeMultiplicity << "\tduplicate\t"
                      << contig_name << '\t' << direction << '\t'
                      << long_edge << '\t' << extension_offset << '\t'
                      << duplicate_source << '\n';
        continue;
      }
      emitted_extensions.Add(key);
      PackKeyAsEdge(key, edge_len, base_meta.words_per_edge,
                    kContigEdgeMultiplicity, &packed);
      writer.Write(packed.data(), extended_edges_added % opt.num_threads,
                   "contig_overlap_extension");
      ++extended_edges_added;
      if (opt.human_readable) contig_report << BasesFromKey(key, edge_len) << '\t'
                    << kContigEdgeMultiplicity << "\tadded\t"
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

  if (!opt.long_contig_file.empty()) {
    std::ifstream long_contigs(opt.long_contig_file);
    std::string name, seq;
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
  writer.Finalize();

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
  summary << "merge_read_count_bloom_bits\t" << read_count_edges.num_bits() << '\n';
  summary << "merge_read_count_bloom_bytes\t" << read_count_edges.num_bytes() << '\n';
  summary << "merge_read_count_bloom_hashes\t" << read_count_edges.num_hashes() << '\n';
  summary << "merge_read_count_bloom_estimated_fpr\t"
          << std::setprecision(12)
          << EstimateFpr(read_count_edges.num_bits(), expected_merge_items,
                         read_count_edges.num_hashes()) << '\n';
  summary << "merge_extension_bloom_bits\t" << emitted_extensions.num_bits() << '\n';
  summary << "merge_extension_bloom_bytes\t" << emitted_extensions.num_bytes() << '\n';
  summary << "merge_extension_bloom_hashes\t" << emitted_extensions.num_hashes() << '\n';
  summary << "merge_extension_bloom_estimated_fpr\t"
          << std::setprecision(12)
          << EstimateFpr(emitted_extensions.num_bits(), expected_merge_items,
                         emitted_extensions.num_hashes()) << '\n';
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
  summary << "contig_edge_multiplicity\t" << kContigEdgeMultiplicity << '\n';
  summary << "output_edges_total\t" << base_added + extended_edges_added << '\n';
    summary << "representation\tbinary_edges_read_count_plus_exact_k_overlap_extensions_bloom_redundancy_filter\n";
  }
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
