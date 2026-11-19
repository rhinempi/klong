#ifndef MEGAHIT_EXACT_MERGE_DEDUP_H
#define MEGAHIT_EXACT_MERGE_DEDUP_H

#include <cstdint>
#include <functional>
#include <memory>
#include <string>

struct ExactDedupStats {
  uint64_t bytes_written, bytes_read, runs_written, buffer_bytes;
};

class ExactMergeDedup {
 public:
  virtual ~ExactMergeDedup() = default;
  virtual void AddBase(const uint64_t *key, int worker) = 0;
  virtual void AddCandidate(const uint64_t *key, uint64_t ordinal) = 0;
  virtual void AddRecord(const uint64_t *key, uint64_t priority, int worker) = 0;
  virtual void FinishUnordered(const std::function<void(const uint64_t *, uint64_t, int)> &emit) = 0;
  virtual void Check() const = 0;
  virtual void Finish(const std::function<void(const uint64_t *, uint64_t)> &emit) = 0;
  virtual ExactDedupStats Stats() const = 0;
};

std::unique_ptr<ExactMergeDedup> MakeExactMergeDedup(unsigned words, int workers,
                                                  uint64_t buffer_bytes,
                                                  const std::string &output_prefix);

#endif
