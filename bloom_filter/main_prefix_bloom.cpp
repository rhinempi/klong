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
#include <thread>
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
constexpr uint64_t kMinReserveBytes = 1ULL * 1024 * 1024 * 1024;
constexpr uint64_t kSplitMixSeed = 0x9e3779b97f4a7c15ULL;

struct Options {
  std::string short_prefix;
  std::string long_prefix;
  std::string output_prefix;
  unsigned short_k{0};
  unsigned long_k{0};
  int num_threads{1};
  double max_bloom_gb{60.0};
  double reserve_gb{1.0};
  int max_hashes{16};
};

struct KmerKey {
  static constexpr int kWords = (kMaxK * 2 + 63) / 64;
  uint64_t words[kWords];

  KmerKey() { std::fill(words, words + kWords, 0); }
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

unsigned EdgeLength(const EdgeIoMetadata &meta) { return meta.kmer_size + 1; }

int CompareForwardAndReverseComplement(const uint32_t *edge, unsigned offset,
                                       unsigned k) {
  for (unsigned i = 0; i < k; ++i) {
    uint8_t f = GetBase(edge, offset + i);
    uint8_t r = 3 - GetBase(edge, offset + k - 1 - i);
    if (f != r) return f < r ? -1 : 1;
  }
  return 0;
}

KmerKey CanonicalKmerFromEdge(const uint32_t *edge, unsigned offset,
                              unsigned k) {
  KmerKey key;
  bool use_forward = CompareForwardAndReverseComplement(edge, offset, k) <= 0;
  for (unsigned i = 0; i < k; ++i) {
    uint8_t base = use_forward ? GetBase(edge, offset + i)
                               : static_cast<uint8_t>(3 - GetBase(edge, offset + k - 1 - i));
    SetPackedBase(&key, i, base);
  }
  return key;
}

std::string BasesFromEdge(const uint32_t *edge, unsigned offset, unsigned k) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  std::string s;
  s.reserve(k);
  for (unsigned i = 0; i < k; ++i) s.push_back(bases[GetBase(edge, offset + i)]);
  return s;
}

std::string ReverseComplementBasesFromEdge(const uint32_t *edge, unsigned offset, unsigned k) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  std::string s;
  s.reserve(k);
  for (unsigned i = 0; i < k; ++i) {
    s.push_back(bases[3 - GetBase(edge, offset + k - 1 - i)]);
  }
  return s;
}

std::string BasesFromKey(const KmerKey &key, unsigned k) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  std::string s;
  s.reserve(k);
  for (unsigned i = 0; i < k; ++i) {
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

EdgeIoMetadata LoadMetadata(const std::string &prefix) {
  std::ifstream is(prefix + ".edges.info");
  if (!is) xfatal("Cannot open {s}.edges.info\n", prefix.c_str());
  EdgeIoMetadata meta;
  meta.Deserialize(is);
  return meta;
}

uint64_t FileRecordCount(const EdgeIoMetadata &meta, unsigned file_id) {
  if (!meta.is_sorted) {
    return file_id == 0 ? static_cast<uint64_t>(meta.num_edges) : 0;
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

void ForEachEdgeInFile(const std::string &prefix, const EdgeIoMetadata &meta,
                       unsigned file_id,
                       const std::function<void(const uint32_t *)> &fn) {
  uint64_t n = FileRecordCount(meta, file_id);
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

void WriteSummary(const std::string &path, const Options &opt,
                  const EdgeIoMetadata &short_meta, const EdgeIoMetadata &long_meta,
                  const BloomFilter &bloom, uint64_t long_records,
                  uint64_t bloom_insertions, uint64_t short_observations,
                  uint64_t bloom_hits, uint64_t bloom_misses, unsigned key_len) {
  std::ofstream out(path);
  out << "short_k\t" << opt.short_k << '\n';
  out << "long_k\t" << opt.long_k << '\n';
  out << "short_edge_prefix\t" << opt.short_prefix << '\n';
  out << "long_edge_prefix\t" << opt.long_prefix << '\n';
  out << "short_edge_kmer_size\t" << short_meta.kmer_size << '\n';
  out << "long_edge_kmer_size\t" << long_meta.kmer_size << '\n';
  out << "long_edge_records\t" << long_records << '\n';
  out << "bloom_insertions\t" << bloom_insertions << '\n';
  out << "reduction_key_bases\t" << key_len << '\n';
  out << "short_edge_observations\t" << short_observations << '\n';
  out << "bloom_hits_probably_removed\t" << bloom_hits << '\n';
  out << "bloom_misses_not_removed\t" << bloom_misses << '\n';
  out << "bloom_bits\t" << bloom.num_bits() << '\n';
  out << "bloom_bytes_each\t" << bloom.num_bytes() << '\n';
  out << "bloom_filters\t2\n";
  out << "bloom_total_bytes\t" << bloom.num_bytes() * 2 << '\n';
  out << "bloom_gb_each\t" << std::setprecision(6)
      << (static_cast<double>(bloom.num_bytes()) / 1024.0 / 1024.0 / 1024.0) << '\n';
  out << "bloom_total_gb\t"
      << (static_cast<double>(bloom.num_bytes()) * 2 / 1024.0 / 1024.0 / 1024.0) << '\n';
  out << "bloom_hashes\t" << bloom.num_hashes() << '\n';
  out << "estimated_false_positive_rate\t"
      << std::setprecision(12) << EstimateFpr(bloom.num_bits(), bloom_insertions, bloom.num_hashes()) << '\n';
  out << "threads\t" << opt.num_threads << '\n';
  out << "representation\tbinary_edges_2bit_canonical_edge_prefixes_both_orientations\n";
  out << "removed_file\t" << opt.output_prefix << ".removed.txt\n";
  out << "non_removed_file\t" << opt.output_prefix << ".non_removed.txt\n";
}

Options ParseOptions(int argc, char **argv) {
  Options opt;
  OptionsDescription desc;
  desc.AddOption("short_prefix", "s", opt.short_prefix, "(*) prefix of shorter k binary edge files");
  desc.AddOption("long_prefix", "l", opt.long_prefix, "(*) prefix of longer k binary edge files");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "(*) output prefix for debug reports");
  desc.AddOption("short_k", "k", opt.short_k, "(*) shorter k-mer size whose (k+1)-edges are compared");
  desc.AddOption("long_k", "K", opt.long_k, "(*) longer k-mer size to report");
  desc.AddOption("num_cpu_threads", "t", opt.num_threads, "number of CPU threads");
  desc.AddOption("max_bloom_gb", "m", opt.max_bloom_gb, "max Bloom filter memory in GB");
  desc.AddOption("reserve_gb", "r", opt.reserve_gb, "RAM to leave unused in GB");
  desc.AddOption("max_hashes", "H", opt.max_hashes, "maximum number of Bloom hash probes");
  desc.Parse(argc, argv);

  if (opt.short_prefix.empty() || opt.long_prefix.empty() || opt.output_prefix.empty() ||
      opt.short_k == 0 || opt.long_k == 0) {
    std::cerr << "Usage: " << argv[0] << " " << std::string(desc) << std::endl;
    std::exit(1);
  }
  if (opt.num_threads <= 0) opt.num_threads = 1;
  if (opt.max_hashes <= 0) opt.max_hashes = 1;
  return opt;
}

}  // namespace

int main_prefix_bloom(int argc, char **argv) {
  Options opt = ParseOptions(argc, argv);
  EdgeIoMetadata short_meta = LoadMetadata(opt.short_prefix);
  EdgeIoMetadata long_meta = LoadMetadata(opt.long_prefix);

  if (short_meta.kmer_size != opt.short_k || long_meta.kmer_size != opt.long_k) {
    xfatal("Invalid k sizes: short edge k {}, long edge k {}, requested short {}, long {}\n",
           short_meta.kmer_size, long_meta.kmer_size, opt.short_k, opt.long_k);
  }

  const unsigned key_len = EdgeLength(short_meta);
  if (EdgeLength(long_meta) < key_len) {
    xfatal("Long edge length {} is shorter than reduction key length {}\n",
           EdgeLength(long_meta), key_len);
  }

  uint64_t long_records = TotalRecords(long_meta);
  uint64_t bloom_insertions = long_records * 2;
  uint64_t bloom_byte_cap = SelectBloomByteCap(opt);
  uint64_t target_bytes = std::max<uint64_t>(64ULL << 20, bloom_insertions * 64ULL);
  uint64_t bloom_bytes = std::min(bloom_byte_cap / 2, target_bytes);
  if (bloom_bytes < (8ULL << 20)) bloom_bytes = std::min<uint64_t>(bloom_byte_cap, 8ULL << 20);
  uint64_t bloom_bits = bloom_bytes * 8;
  int hashes = ChooseHashCount(bloom_bits, std::max<uint64_t>(1, bloom_insertions), opt.max_hashes);
  BloomFilter bloom(bloom_bits, hashes);
  BloomFilter candidate_bloom(bloom_bits, hashes);

  xinfo("Prefix Bloom: long records {}, Bloom bytes each {}, hashes {}, estimated FPR {}\n",
        long_records, bloom.num_bytes(), bloom.num_hashes(),
        EstimateFpr(bloom.num_bits(), bloom_insertions, bloom.num_hashes()));

#pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(long_meta.num_files); ++file_id) {
    ForEachEdgeInFile(opt.long_prefix, long_meta, file_id, [&](const uint32_t *edge) {
      bloom.Add(CanonicalKmerFromEdge(edge, 0, key_len));
      bloom.Add(CanonicalKmerFromEdge(edge, EdgeLength(long_meta) - key_len, key_len));
    });
  }

  std::ofstream non_removed(opt.output_prefix + ".non_removed.txt");
  non_removed << "#canonical_short_edge\tshort_edge_as_seen\tshort_edge_multiplicity\treason\n";

  std::atomic<uint64_t> short_observations{0};
  std::atomic<uint64_t> bloom_hits{0};
  std::atomic<uint64_t> bloom_misses{0};
  std::vector<std::ostringstream> non_removed_buffers(opt.num_threads);

#pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(short_meta.num_files); ++file_id) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    ForEachEdgeInFile(opt.short_prefix, short_meta, file_id, [&](const uint32_t *edge) {
      KmerKey key = CanonicalKmerFromEdge(edge, 0, key_len);
      bool hit = bloom.Contains(key);
      short_observations.fetch_add(1, std::memory_order_relaxed);
      if (hit) {
        candidate_bloom.Add(key);
        bloom_hits.fetch_add(1, std::memory_order_relaxed);
      } else {
        non_removed_buffers[tid] << BasesFromKey(key, key_len) << '\t'
                                 << BasesFromEdge(edge, 0, key_len) << '\t'
                                 << Multiplicity(edge, short_meta.words_per_edge) << '\t'
                                 << "bloom_absent" << '\n';
        bloom_misses.fetch_add(1, std::memory_order_relaxed);
      }
    });
  }

  for (int i = 0; i < opt.num_threads; ++i) {
    non_removed << non_removed_buffers[i].str();
  }

  std::ofstream removed(opt.output_prefix + ".removed.txt");
  removed << "#canonical_short_edge\tlong_edge_prefix\tcanonical_long_edge_prefix\tlong_edge\tlong_edge_multiplicity\torientation\treason\n";
  std::vector<std::ostringstream> removed_buffers(opt.num_threads);

#pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(long_meta.num_files); ++file_id) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    ForEachEdgeInFile(opt.long_prefix, long_meta, file_id, [&](const uint32_t *edge) {
      auto emit_long_prefix = [&](unsigned offset, const char *orientation) {
        KmerKey short_key = CanonicalKmerFromEdge(edge, offset, key_len);
        if (!candidate_bloom.Contains(short_key)) return;
        std::string oriented_prefix =
            offset == 0 ? BasesFromEdge(edge, offset, key_len)
                        : ReverseComplementBasesFromEdge(edge, offset, key_len);
        removed_buffers[tid] << BasesFromKey(short_key, key_len) << '\t'
                             << oriented_prefix << '\t'
                             << BasesFromKey(short_key, key_len) << '\t'
                             << BasesFromEdge(edge, 0, EdgeLength(long_meta)) << '\t'
                             << Multiplicity(edge, long_meta.words_per_edge) << '\t'
                             << orientation << '\t'
                             << "candidate_bloom_present_long_edge_prefix_counterpart" << '\n';
      };
      emit_long_prefix(0, "forward_prefix");
      emit_long_prefix(EdgeLength(long_meta) - key_len, "reverse_complement_prefix");
    });
  }
  for (int i = 0; i < opt.num_threads; ++i) removed << removed_buffers[i].str();

  WriteSummary(opt.output_prefix + ".summary.txt", opt, short_meta, long_meta, bloom,
               long_records, bloom_insertions, short_observations.load(),
               bloom_hits.load(), bloom_misses.load(), key_len);

  xinfo("Prefix Bloom report wrote {s}.removed.txt and {s}.non_removed.txt\n",
        opt.output_prefix.c_str(), opt.output_prefix.c_str());
  return 0;
}
