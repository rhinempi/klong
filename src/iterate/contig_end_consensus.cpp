#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "sequence/io/fastx_reader.h"
#include "sequence/sequence_package.h"

namespace {
using Counts = std::array<uint64_t, 4>;

struct Reference {
  std::string name;
  std::string sequence;
  std::vector<Counts> left, right;
  uint64_t primary_reads{0};
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

void ProcessSam(const std::string &path, std::vector<Reference> *references,
                const std::unordered_map<std::string, size_t> &index) {
  std::ifstream source(path);
  if (!source) throw std::runtime_error("Cannot open SAM: " + path);
  SeqPackage packed_read;
  std::string line;
  uint64_t line_number = 0;
  while (std::getline(source, line)) {
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
      uint64_t flag = Number(fields[1]);
      if (flag & (4 | 256 | 2048)) continue;
      auto entry = index.find(fields[2]);
      if (entry == index.end()) throw std::runtime_error("Unknown SAM reference: " + fields[2]);
      auto &ref = (*references)[entry->second];
      ++ref.primary_reads;
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
      ref.overhang_reads += voted;
    } catch (const std::exception &error) {
      throw std::runtime_error(path + ":" + std::to_string(line_number) + ": " + error.what());
    }
  }
  if (!source.eof()) throw std::runtime_error("Error reading SAM: " + path);
}
}  // namespace

int main_contig_end_consensus(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "Usage: endconsensus <ends.fa> <extended.fa> <stats.tsv> <alignments.sam>...\n";
    return 1;
  }
  try {
    std::vector<Reference> references;
    std::unordered_map<std::string, size_t> index;
    FastxReader reader(argv[1]);
    while (auto record = reader.ReadNext()) {
      Reference ref;
      ref.name.assign(record->name.s, record->name.l);
      ref.sequence.assign(record->seq.s, record->seq.l);
      if (!index.emplace(ref.name, references.size()).second) throw std::runtime_error("Duplicate reference: " + ref.name);
      ref.extend_left = !EndsWith(ref.name, "|end");
      ref.extend_right = !EndsWith(ref.name, "|start");
      references.emplace_back(std::move(ref));
    }
    for (int i = 4; i < argc; ++i) ProcessSam(argv[i], &references, index);
    std::ofstream output(argv[2]), stats(argv[3]), info(std::string(argv[2]) + ".info");
    if (!output || !stats || !info) throw std::runtime_error("Cannot open consensus output files");
    stats << "reference\treference_length\tleft_added\tright_added\tprimary_reads\toverhang_reads\n";
    uint64_t extended = 0, added = 0, output_bases = 0;
    for (const auto &ref : references) {
      auto left = Consensus(ref.left), right = Consensus(ref.right);
      std::reverse(left.begin(), left.end());
      stats << ref.name << '\t' << ref.sequence.size() << '\t' << left.size() << '\t'
            << right.size() << '\t' << ref.primary_reads << '\t' << ref.overhang_reads << '\n';
      if (left.empty() && right.empty()) continue;
      auto length = left.size() + ref.sequence.size() + right.size();
      output << '>' << ref.name << " flag=0 multi=1 reference_length=" << ref.sequence.size()
             << " left_extension=" << left.size() << " right_extension=" << right.size()
             << " len=" << length << '\n' << left << ref.sequence << right << '\n';
      ++extended;
      added += left.size() + right.size();
      output_bases += length;
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
