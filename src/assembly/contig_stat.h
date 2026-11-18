//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_CONTIG_STAT_H
#define MEGAHIT_CONTIG_STAT_H

#include <map>
#include <omp.h>
#include <vector>
#include "unitig_graph.h"
#include "utils/histgram.h"
#include "utils/utils.h"

using ContigStat = std::map<std::string, uint64_t>;

inline ContigStat CalcAndPrintStat(UnitigGraph &graph, bool print = true,
                                   bool changed_only = false) {
  uint32_t n_isolated = 0, n_looped = 0;
  Histgram<uint64_t> hist;
  std::vector<std::map<uint64_t, size_t>> local_counts(omp_get_max_threads());
#pragma omp parallel reduction(+ : n_looped, n_isolated)
  {
    std::map<uint64_t, size_t> counts;
#pragma omp for
    for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
      auto adapter = graph.MakeVertexAdapter(i);
      if (changed_only && !adapter.IsChanged()) {
        continue;
      }
      ++counts[adapter.GetLength() + graph.k()];
      n_looped += adapter.IsLoop();
      n_isolated += adapter.IsStandalone() || (graph.InDegree(adapter) == 0 &&
                                               graph.OutDegree(adapter) == 0);
    }
    local_counts[omp_get_thread_num()].swap(counts);
  }
  for (const auto &counts : local_counts)
    for (const auto &entry : counts) hist.insert(entry.first, entry.second);
  uint64_t total_size = hist.sum();
  ContigStat stat = {{"Max", hist.maximum()},
                     {"Min", hist.minimum()},
                     {"N50", hist.Nx(0.5 * total_size)},
                     {"total size", total_size},
                     {"number contigs", hist.size()},
                     {"number looped", n_looped},
                     {"number isolated", n_isolated}};

  if (print) {
    xinfo("");
    for (auto &kv : stat) {
      xinfoc("{s}: {}, ", kv.first.c_str(), kv.second);
    }
    xinfoc("{s}", "\n");
  }
  return stat;
}

#endif  // MEGAHIT_CONTIG_STAT_H
