#ifndef MEGAHIT_ANCHOR_PATHS_H
#define MEGAHIT_ANCHOR_PATHS_H

#include <cstddef>
#include <string>

class SDBG;

constexpr size_t kFlankingAnchorEndLength = 290;

struct AnchorPathStats {
  size_t queries{0}, found{0}, crossing_branches{0};
  size_t extended{0}, added_nt{0};
};

AnchorPathStats WriteAnchorPathReport(const SDBG &graph,
                                      const std::string &source,
                                      size_t end_length,
                                      const std::string &output,
                                      const std::string &threaded_fasta,
                                      bool stitched_only);

#endif  // MEGAHIT_ANCHOR_PATHS_H
