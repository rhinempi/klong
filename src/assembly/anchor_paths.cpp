#include "anchor_paths.h"

#include <omp.h>
#include <algorithm>
#include <cstdio>
#include <exception>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "definitions.h"
#include "sdbg/sdbg.h"
#include "sequence/io/fastx_reader.h"

namespace {

enum class Match { Full, Missing, TooShort, Circular, NonAcgt, RawUnitig };

const char *MatchName(Match match) {
  switch (match) {
    case Match::Full: return "full";
    case Match::Missing: return "missing";
    case Match::TooShort: return "too_short";
    case Match::Circular: return "circular";
    case Match::NonAcgt: return "non_acgt";
    case Match::RawUnitig: return "raw_unitig";
  }
  return "unknown";
}

struct Query {
  std::string name, anchor;
  Match match{Match::Missing};
  size_t matched_nt{0}, branches{0};
  std::string threaded;
  size_t left_nt{0}, right_nt{0};
  std::exception_ptr error;
};

bool IsCircular(const char *comment) {
  std::istringstream fields(comment ? comment : "");
  std::string token;
  while (fields >> token) {
    if (token.compare(0, 5, "flag=") == 0)
      return (std::stoul(token.substr(5)) & contig_flag::kLoop) != 0;
  }
  return false;
}

uint8_t BaseCode(char base) {
  switch (base) {
    case 'A': case 'a': return 1;
    case 'C': case 'c': return 2;
    case 'G': case 'g': return 3;
    case 'T': case 't': return 4;
    default: return 0;
  }
}

uint8_t EdgeBase(const SDBG &graph, uint64_t edge) {
  auto base = graph.GetW(edge);
  return base > kAlphabetSize ? base - kAlphabetSize : base;
}

void SpellPath(const SDBG &graph, const std::vector<uint64_t> &path, Query &query) {
  constexpr char dna[] = "NACGT";
  std::unordered_set<uint64_t> used(path.begin(), path.end());
  std::string left_steps, right_steps;
  auto first = path.front();
  auto prev = graph.PrevSimplePathEdge(first);
  while (prev != SDBG::kNullID && used.insert(prev).second) {
    first = prev;
    left_steps.push_back(dna[EdgeBase(graph, first)]);
    prev = graph.PrevSimplePathEdge(first);
  }
  auto next = graph.NextSimplePathEdge(path.back());
  while (next != SDBG::kNullID && used.insert(next).second) {
    right_steps.push_back(dna[EdgeBase(graph, next)]);
    next = graph.NextSimplePathEdge(next);
  }

  query.left_nt = left_steps.size();
  query.right_nt = right_steps.size();
  query.threaded.reserve(query.anchor.size() + query.left_nt + query.right_nt);
  uint8_t seed[kMaxK];
  graph.GetLabel(first, seed);
  for (size_t i = 0; i < graph.k(); ++i) query.threaded.push_back(dna[seed[i]]);
  query.threaded.append(left_steps.rbegin(), left_steps.rend());
  for (auto edge : path) query.threaded.push_back(dna[EdgeBase(graph, edge)]);
  query.threaded += right_steps;
}

void FindPath(const SDBG &graph, Query &query, bool emit_contigs) {
  if (query.match != Match::Missing) return;
  const auto k = graph.k();
  if (query.anchor.size() <= k) {
    query.match = Match::TooShort;
    return;
  }
  uint8_t bases[2 * kFlankingAnchorEndLength];
  for (size_t i = 0; i < query.anchor.size(); ++i) {
    bases[i] = BaseCode(query.anchor[i]);
    if (!bases[i]) {
      query.match = Match::NonAcgt;
      return;
    }
  }
  if (!graph.size()) return;
  auto node_end = graph.IndexBinarySearch(bases);
  std::vector<uint64_t> path;
  if (emit_contigs && node_end != SDBG::kNullID) path.reserve(query.anchor.size() - k);
  for (size_t i = k; i < query.anchor.size() && node_end != SDBG::kNullID; ++i) {
    // Walk the outgoing edge group even if its last slot was pruned.
    uint64_t next = SDBG::kNullID;
    unsigned outdegree = 0;
    auto edge = node_end;
    do {
      if (graph.IsValidEdge(edge)) {
        ++outdegree;
        if (EdgeBase(graph, edge) == bases[i]) next = edge;
      }
      --edge;
    } while (edge != SDBG::kNullID && !graph.IsLastOrTip(edge));
    if (next == SDBG::kNullID) return;
    if (i > k && (outdegree > 1 || graph.EdgeIndegree(next) > 1)) ++query.branches;
    query.matched_nt = i + 1;
    if (emit_contigs) path.push_back(next);
    if (i + 1 == query.anchor.size()) {
      query.match = Match::Full;
      if (emit_contigs) SpellPath(graph, path, query);
      return;
    }
    node_end = graph.Forward(next);
  }
}

}  // namespace

AnchorPathStats WriteAnchorPathReport(const SDBG &graph,
                                      const std::string &source,
                                      size_t end_length,
                                      const std::string &output,
                                      const std::string &threaded_fasta,
                                      bool stitched_only) {
  if (end_length < graph.k() || end_length > kFlankingAnchorEndLength)
    throw std::invalid_argument("Invalid end length for anchor path report");
  FastxReader reader(source);
  const auto temp = output + ".tmp";
  std::ofstream report(temp);
  if (!report) throw std::runtime_error("Cannot open anchor path report: " + temp);
  report << "source_contig\tanchor_nt\tmatch\tmatched_prefix_nt\tbranch_crossings\n";
  const bool emit_contigs = !threaded_fasta.empty();
  const auto fasta_temp = threaded_fasta + ".tmp";
  std::ofstream contigs;
  if (emit_contigs) {
    contigs.open(fasta_temp);
    if (!contigs) throw std::runtime_error("Cannot open threaded contig FASTA: " + fasta_temp);
  }

  AnchorPathStats stats;
  std::vector<Query> batch;
  const size_t batch_size = emit_contigs ? 256 : 4096;
  batch.reserve(batch_size);
  const auto flush = [&]() {
    if (batch.empty()) return;
    const int workers = std::min(omp_get_max_threads(), static_cast<int>(batch.size()));
#pragma omp parallel for num_threads(workers) if(workers > 1)
    for (size_t i = 0; i < batch.size(); ++i) {
      try {
        FindPath(graph, batch[i], emit_contigs);
      } catch (...) {
        batch[i].error = std::current_exception();
      }
    }
    for (const auto &query : batch) {
      if (query.error) std::rethrow_exception(query.error);
      if (query.match == Match::Full || query.match == Match::Missing) ++stats.queries;
      if (query.match == Match::Full) {
        ++stats.found;
        if (query.branches) ++stats.crossing_branches;
        if (emit_contigs) {
          const auto added = query.left_nt + query.right_nt;
          stats.added_nt += added;
          stats.extended += added != 0;
          contigs << '>' << query.name << " anchor_nt=" << query.anchor.size()
                  << " left_extension=" << query.left_nt << " right_extension=" << query.right_nt
                  << " branch_crossings=" << query.branches << " len=" << query.threaded.size()
                  << '\n' << query.threaded << '\n';
        }
      }
      report << query.name << '\t' << query.anchor.size() << '\t' << MatchName(query.match)
             << '\t' << query.matched_nt << '\t' << query.branches << '\n';
    }
    batch.clear();
  };
  while (auto record = reader.ReadNext()) {
    batch.emplace_back();
    auto &query = batch.back();
    query.name.assign(record->name.s, record->name.l);
    if (IsCircular(record->comment.l ? record->comment.s : nullptr)) {
      query.match = Match::Circular;
    } else if (record->seq.l <= 2 * end_length) {
      query.anchor.assign(record->seq.s, record->seq.l);
      if (stitched_only) query.match = Match::RawUnitig;
    } else {
      query.anchor.assign(record->seq.s, end_length);
      query.anchor.append(record->seq.s + record->seq.l - end_length, end_length);
    }
    if (batch.size() == batch_size) flush();
  }
  flush();
  if (emit_contigs) {
    contigs.close();
    if (!contigs || std::rename(fasta_temp.c_str(), threaded_fasta.c_str()) != 0)
      throw std::runtime_error("Cannot write threaded contig FASTA: " + threaded_fasta);
  }
  report.close();
  if (!report || std::rename(temp.c_str(), output.c_str()) != 0)
    throw std::runtime_error("Cannot write anchor path report: " + output);
  return stats;
}
