#ifndef MEGAHIT_NATIVE_END_ANCHORS_H
#define MEGAHIT_NATIVE_END_ANCHORS_H

#include <cstddef>
#include <string>

namespace native {
void write_end_anchors(const std::string &source, const std::string &output,
                       size_t end_length);
void write_outward_consensus(const std::string &source, const std::string &output,
                             size_t end_length);

struct AnchorRestoreStats {
  size_t contigs{0}, joined{0}, extended{0}, bases_added{0};
};
AnchorRestoreStats restore_end_anchors(const std::string &source, const std::string &unitigs,
                                       const std::string &output, size_t end_length, size_t k, int threads,
                                       int max_mismatches = -1, const std::string &report = "",
                                       const std::string &short_unitigs = "", bool best_hit = false,
                                       int min_identity = 0);
}  // namespace native

#endif
