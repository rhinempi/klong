#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "sequence/io/fastx_reader.h"
#include "sequence/sequence_package.h"
#include "sdbg/sdbg_def.h"

namespace {
using Counts = std::array<uint64_t, 4>;
using VoteLocks = std::array<std::mutex, 4096>;

struct SamRange {
  std::string path;
  std::streamoff begin, end;
};

struct Reference {
  std::string name;
  std::string sequence;
  std::vector<Counts> left, right;
  uint64_t primary_reads{0};
  uint64_t secondary_reads{0};
  uint64_t overhang_reads{0};
  bool extend_left{true};
  bool extend_right{true};
};

struct CigarOp {
  uint32_t length;
  char code;
};

bool EndsWith(const std::string &value, const std::string &suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> Fields(const std::string &line) {
  std::vector<std::string> result;
  size_t begin = 0;
  while (true) {
    auto end = line.find('\t', begin);
    result.emplace_back(line.substr(begin, end - begin));
    if (end == std::string::npos) return result;
    begin = end + 1;
  }
}

uint64_t Number(const std::string &value) {
  if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos) {
    throw std::runtime_error("Invalid SAM integer: " + value);
  }
  return std::stoull(value);
}

std::vector<CigarOp> ParseCigar(const std::string &text) {
  std::vector<CigarOp> result;
  uint64_t length = 0;
  for (char c : text) {
    if (c >= '0' && c <= '9') {
      length = length * 10 + c - '0';
      if (length > std::numeric_limits<int32_t>::max()) {
        throw std::runtime_error("CIGAR operation is too long");
      }
    } else {
      if (length == 0 || std::string("MIDNSHP=X").find(c) == std::string::npos) {
        throw std::runtime_error("Invalid CIGAR: " + text);
      }
      result.push_back({static_cast<uint32_t>(length), c});
      length = 0;
    }
  }
  if (length || result.empty() || (text.back() >= '0' && text.back() <= '9')) {
    throw std::runtime_error("Invalid CIGAR: " + text);
  }
  return result;
}

bool KnownBase(char c) {
  return c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
         c == 'a' || c == 'c' || c == 'g' || c == 't';
}

void AddVote(std::vector<Counts> *votes, size_t distance, uint8_t base) {
  if (votes->size() <= distance) votes->resize(distance + 1, Counts{});
  ++(*votes)[distance][base];
}

std::string Consensus(const std::vector<Counts> &votes) {
  std::string result;
  for (const auto &counts : votes) {
    auto best = std::max_element(counts.begin(), counts.end());
    if (*best == 0 || std::count(counts.begin(), counts.end(), *best) != 1) break;
    result.push_back("ACGT"[best - counts.begin()]);
  }
  return result;
}

// Move each cut past its current query group, preserving per-read deduplication.
std::streamoff QueryBoundary(const std::string &path, std::streamoff offset,
                             std::streamoff length) {
  if (offset == 0 || offset == length) return offset;
  std::ifstream source(path, std::ios::binary);
  if (!source) throw std::runtime_error("Cannot open SAM: " + path);
  source.seekg(offset - 1);
  char previous;
  if (!source.get(previous)) throw std::runtime_error("Cannot seek SAM: " + path);
  std::string line, query;
  if (previous != '\n') std::getline(source, line);
  bool first = true;
  while (source) {
    auto position = source.tellg();
    if (!std::getline(source, line)) break;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line[0] == '@') return position;
    auto name = line.substr(0, line.find('\t'));
    if (!first && name != query) return position;
    query = std::move(name);
    first = false;
  }
  if (source.bad()) throw std::runtime_error("Error reading SAM: " + path);
  return length;
}

std::vector<SamRange> SplitSam(const std::vector<std::string> &paths, int threads) {
  std::vector<SamRange> ranges;
  for (const auto &path : paths) {
    std::ifstream source(path, std::ios::binary | std::ios::ate);
    if (!source) throw std::runtime_error("Cannot open SAM: " + path);
    std::streamoff length = source.tellg();
    if (length < 0) throw std::runtime_error("Cannot determine SAM size: " + path);
    std::streamoff begin = 0;
    for (int i = 1; i <= threads; ++i) {
      auto target = (length / threads) * i + (length % threads) * i / threads;
      auto end = QueryBoundary(path, target, length);
      if (end > begin) ranges.push_back({path, begin, end});
      begin = end;
    }
  }
  return ranges;
}

void ProcessSam(const SamRange &range, std::vector<Reference> *references,
                const std::unordered_map<std::string, size_t> &index, VoteLocks &locks) {
  const auto &path = range.path;
  std::ifstream source(path, std::ios::binary);
  if (!source) throw std::runtime_error("Cannot open SAM: " + path);
  source.seekg(range.begin);
  if (!source) throw std::runtime_error("Cannot seek SAM: " + path);
  SeqPackage packed_read;
  std::string current_read;
  std::vector<size_t> voted_references;
  std::string line;
  uint64_t line_number = 0;
  std::streamoff position_in_file = range.begin;
  while (position_in_file < range.end && std::getline(source, line)) {
    auto line_offset = position_in_file;
    position_in_file += static_cast<std::streamoff>(line.size()) + (source.eof() ? 0 : 1);
    ++line_number;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    try {
      auto fields = Fields(line);
      if (line[0] == '@') {
        if (fields[0] == "@SQ") {
          std::string name, length;
          for (const auto &field : fields) {
            if (field.compare(0, 3, "SN:") == 0) name = field.substr(3);
            if (field.compare(0, 3, "LN:") == 0) length = field.substr(3);
          }
          auto entry = index.find(name);
          if (entry == index.end() || Number(length) != (*references)[entry->second].sequence.size()) {
            throw std::runtime_error("SAM reference header does not match end FASTA: " + name);
          }
        }
        continue;
      }
      if (fields.size() < 11) throw std::runtime_error("SAM requires at least 11 fields");
      // Minimap2 groups a query's alignments; retain only that query's voted ends.
      if (fields[0] != current_read) {
        current_read = fields[0];
        voted_references.clear();
      }
      uint64_t flag = Number(fields[1]);
      if (flag & (4 | 256 | 2048)) continue;
      auto entry = index.find(fields[2]);
      if (entry == index.end()) throw std::runtime_error("Unknown SAM reference: " + fields[2]);
      auto &ref = (*references)[entry->second];
      auto &mutex = locks[entry->second % locks.size()];
      {
        std::lock_guard<std::mutex> lock(mutex);
        ++ref.primary_reads;
      }
      if (std::find(voted_references.begin(), voted_references.end(), entry->second) != voted_references.end()) continue;
      if (fields[9] == "*") continue;
      uint64_t position = Number(fields[3]);
      if (position == 0 || position > ref.sequence.size()) throw std::runtime_error("SAM POS outside reference");
      auto cigar = ParseCigar(fields[5]);
      uint64_t query_length = 0, reference_length = 0, aligned_length = 0;
      size_t first = 0, last = cigar.size() - 1;
      while (first < cigar.size() && cigar[first].code == 'H') ++first;
      while (last > first && cigar[last].code == 'H') --last;
      for (size_t i = 0; i < cigar.size(); ++i) {
        const auto &op = cigar[i];
        if (op.code == 'M' || op.code == '=' || op.code == 'X' || op.code == 'I' || op.code == 'S') {
          query_length += op.length;
        }
        if (op.code == 'M' || op.code == '=' || op.code == 'X' || op.code == 'D' || op.code == 'N') {
          reference_length += op.length;
        }
        if (op.code == 'M' || op.code == '=' || op.code == 'X') aligned_length += op.length;
        if (op.code == 'S' && i != first && i != last) throw std::runtime_error("Internal soft clip");
        if (op.code == 'H' && i >= first && i <= last) throw std::runtime_error("Internal hard clip");
      }
      if (query_length != fields[9].size() || query_length > std::numeric_limits<unsigned>::max() ||
          aligned_length == 0 || reference_length > ref.sequence.size() - (position - 1)) {
        throw std::runtime_error("CIGAR length does not match SEQ or reference bounds");
      }
      if (fields[10] != "*" && fields[10].size() != query_length) throw std::runtime_error("SEQ/QUAL length mismatch");
      bool has_clip = false;
      for (const auto &op : cigar) has_clip |= op.code == 'S';
      if (!has_clip) continue;

      // SAM SEQ is already oriented to the reference, including FLAG 16 reads.
      packed_read.Clear();
      packed_read.AppendStringSequence(fields[9].data(), static_cast<unsigned>(query_length));
      auto view = packed_read.GetSeqView(0);
      size_t qpos = 0;
      int64_t rpos = static_cast<int64_t>(position - 1);
      bool voted = false;
      std::lock_guard<std::mutex> vote_lock(mutex);
      for (size_t i = 0; i < cigar.size(); ++i) {
        const auto &op = cigar[i];
        if (op.code == 'S') {
          int64_t start = i == first ? rpos - op.length : rpos;
          for (uint32_t j = 0; j < op.length; ++j) {
            int64_t coordinate = start + j;
            if (!KnownBase(fields[9][qpos + j])) continue;
            uint8_t base = view.base_at(qpos + j);
            if (coordinate < 0 && ref.extend_left) {
              AddVote(&ref.left, static_cast<size_t>(-coordinate - 1), base);
              voted = true;
            } else if (coordinate >= static_cast<int64_t>(ref.sequence.size()) && ref.extend_right) {
              AddVote(&ref.right, static_cast<size_t>(coordinate - ref.sequence.size()), base);
              voted = true;
            }
          }
          qpos += op.length;
        } else if (op.code == 'M' || op.code == '=' || op.code == 'X') {
          qpos += op.length;
          rpos += op.length;
        } else if (op.code == 'I') {
          qpos += op.length;
        } else if (op.code == 'D' || op.code == 'N') {
          rpos += op.length;
        }
      }
      if (voted) {
        ++ref.overhang_reads;
        voted_references.push_back(entry->second);
      }
    } catch (const std::exception &error) {
      auto location = range.begin == 0 ? std::to_string(line_number) : "byte " + std::to_string(line_offset);
      throw std::runtime_error(path + ":" + location + ": " + error.what());
    }
  }
  if (source.bad() || (!source.eof() && position_in_file < range.end))
    throw std::runtime_error("Error reading SAM: " + path);
}
}  // namespace

int main_contig_end_consensus(int argc, char **argv) {
  try {
    int first = 1, threads = 1;
    while (first < argc) {
      std::string option(argv[first]);
      if (option == "--") { ++first; break; }
      if (option != "-t" && option != "--threads" && option != "--num-cpu-threads") break;
      if (++first >= argc) throw std::runtime_error("Missing consensus thread count");
      auto value = Number(argv[first++]);
      if (value == 0 || value > static_cast<uint64_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Consensus thread count must be positive and fit an integer");
      threads = static_cast<int>(value);
    }
    if (argc - first < 4) throw std::runtime_error(
        "Usage: endconsensus [-t threads] <ends.fa> <extended.fa> <stats.tsv> <alignments.sam>...");
    threads = std::min(threads, std::max(1, omp_get_num_procs()));
    std::vector<Reference> references;
    std::unordered_map<std::string, size_t> index;
    FastxReader reader(argv[first]);
    while (auto record = reader.ReadNext()) {
      Reference ref;
      ref.name.assign(record->name.s, record->name.l);
      ref.sequence.assign(record->seq.s, record->seq.l);
      if (!index.emplace(ref.name, references.size()).second) throw std::runtime_error("Duplicate reference: " + ref.name);
      ref.extend_left = !EndsWith(ref.name, "|end");
      ref.extend_right = !EndsWith(ref.name, "|start");
      references.emplace_back(std::move(ref));
    }
    std::vector<std::string> paths;
    for (int i = first + 3; i < argc; ++i) paths.emplace_back(argv[i]);
    auto ranges = SplitSam(paths, threads);
    int workers = std::min(threads, std::max(1, static_cast<int>(ranges.size())));
    std::cerr << "Consensus workers: " << workers << ", SAM ranges: " << ranges.size() << '\n';
    VoteLocks locks;
    std::vector<std::exception_ptr> errors(ranges.size());
#pragma omp parallel for num_threads(workers) schedule(dynamic, 1) if(workers > 1)
    for (size_t i = 0; i < ranges.size(); ++i) {
      try { ProcessSam(ranges[i], &references, index, locks); }
      catch (...) { errors[i] = std::current_exception(); }
    }
    for (const auto &error : errors) if (error) std::rethrow_exception(error);
    std::ofstream output(argv[first + 1]), stats(argv[first + 2]), info(std::string(argv[first + 1]) + ".info");
    if (!output || !stats || !info) throw std::runtime_error("Cannot open consensus output files");
    stats << "reference\treference_length\tleft_added\tright_added\tprimary_reads\toverhang_reads\tsecondary_reads\n";
    uint64_t extended = 0, added = 0, output_bases = 0;
    struct Extension { std::string left, right; };
    const size_t block_size = 1024;
    for (size_t begin = 0; begin < references.size(); begin += block_size) {
      size_t count = std::min(block_size, references.size() - begin);
      std::vector<Extension> extensions(count);
      std::vector<std::exception_ptr> failures(count);
#pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1)
      for (size_t i = 0; i < count; ++i) {
        try {
          extensions[i].left = Consensus(references[begin + i].left);
          extensions[i].right = Consensus(references[begin + i].right);
          std::reverse(extensions[i].left.begin(), extensions[i].left.end());
        } catch (...) { failures[i] = std::current_exception(); }
      }
      for (const auto &error : failures) if (error) std::rethrow_exception(error);
      for (size_t i = 0; i < count; ++i) {
        const auto &ref = references[begin + i];
        const auto &left = extensions[i].left;
        const auto &right = extensions[i].right;
        stats << ref.name << '\t' << ref.sequence.size() << '\t' << left.size() << '\t'
              << right.size() << '\t' << ref.primary_reads << '\t' << ref.overhang_reads
              << '\t' << ref.secondary_reads << '\n';
        if (left.empty() && right.empty()) continue;
        auto length = left.size() + ref.sequence.size() + right.size();
        auto multiplicity = std::min(ref.overhang_reads, static_cast<uint64_t>(kMaxMul));
        output << '>' << ref.name << " flag=0 multi=" << multiplicity
               << " support_reads=" << ref.overhang_reads << " reference_length=" << ref.sequence.size()
               << " left_extension=" << left.size() << " right_extension=" << right.size()
               << " len=" << length << '\n' << left << ref.sequence << right << '\n';
        ++extended;
        added += left.size() + right.size();
        output_bases += length;
      }
    }
    info << extended << ' ' << output_bases << '\n';
    output.close();
    stats.close();
    info.close();
    if (!output || !stats || !info) throw std::runtime_error("Error writing consensus output");
    std::cerr << "Consensus: " << references.size() << " references, " << extended
              << " extended, " << added << " bases added\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "endconsensus: " << error.what() << '\n';
    return 1;
  }
}
