#include "end_fm_index.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include "libsais.h"
#include "libsais64.h"

void EndFmIndex::RankBits::Init(size_t length) {
  bits.assign((length + 63) / 64, 0);
}

void EndFmIndex::RankBits::Build() {
  counts.resize(bits.size() + 1);
  counts[0] = 0;
  for (size_t i = 0; i < bits.size(); ++i)
    counts[i + 1] = counts[i] + __builtin_popcountll(bits[i]);
}

uint64_t EndFmIndex::RankBits::Rank(size_t end) const {
  size_t word = end / 64, offset = end % 64;
  return counts[word] + (offset ? __builtin_popcountll(bits[word] & ((uint64_t{1} << offset) - 1)) : 0);
}

void EndFmIndex::Add(const uint64_t *words) {
  for (unsigned i = 0; i < length_; ++i)
    text_.push_back(2 + ((words[i / 32] >> (62 - 2 * (i % 32))) & 3));
  text_.push_back(1);
}

void EndFmIndex::Build(int threads) {
  if (text_.empty()) return;
  // Separators stop cross-contig matches; the unique final sentinel makes LF
  // mapping on the suffix-array-derived BWT unambiguous.
  text_.push_back(0);
  text_size_ = text_.size();
  if (text_size_ <= static_cast<size_t>(std::numeric_limits<int32_t>::max())) {
    std::vector<int32_t> sa(text_size_);
    if (libsais_omp(text_.data(), sa.data(), static_cast<int32_t>(text_size_), 0, nullptr, threads) != 0)
      throw std::runtime_error("libsais failed to build contig-end suffix array");
    Finish(sa, threads);
  } else {
    if (text_size_ > static_cast<size_t>(std::numeric_limits<int64_t>::max()))
      throw std::runtime_error("Contig-end FM index exceeds 64-bit suffix-array capacity");
    std::vector<int64_t> sa(text_size_);
    if (libsais64_omp(text_.data(), sa.data(), static_cast<int64_t>(text_size_), 0, nullptr, threads) != 0)
      throw std::runtime_error("libsais64 failed to build contig-end suffix array");
    Finish(sa, threads);
  }
  std::vector<uint8_t>().swap(text_);
  BuildSeeds(threads);
}

template <typename Index>
void EndFmIndex::Finish(const std::vector<Index> &sa, int threads) {
  for (auto &rank : bwt_) rank.Init(text_size_);
  key_starts_.Init(text_size_);
  size_t previous = Missing();
  for (size_t row = 0; row < text_size_; ++row) {
    size_t position = sa[row];
    if (position == text_size_ - 1 || position % (length_ + 1) != 0) continue;
    // Equal fixed-length records share one query ID, even when repeated in
    // the contig input. Query IDs follow suffix-array (lexicographic) order.
    if (previous == Missing() || std::memcmp(text_.data() + previous, text_.data() + position, length_) != 0) {
      key_starts_.bits[row / 64] |= uint64_t{1} << (row % 64);
      ++keys_;
    }
    previous = position;
  }
#pragma omp parallel for schedule(static) num_threads(threads) if(text_size_ >= 65536)
  for (size_t word = 0; word < bwt_[0].bits.size(); ++word) {
    uint64_t masks[4]{};
    size_t end = std::min(text_size_, word * 64 + 64);
    for (size_t row = word * 64; row < end; ++row) {
      size_t suffix = sa[row];
      uint8_t base = text_[suffix ? suffix - 1 : text_size_ - 1];
      if (base >= 2) masks[base - 2] |= uint64_t{1} << (row % 64);
    }
    for (unsigned base = 0; base < 4; ++base) bwt_[base].bits[word] = masks[base];
  }
  uint64_t start = 1 + (text_size_ - 1) / (length_ + 1);
  for (unsigned base = 0; base < 4; ++base) {
    bwt_[base].Build();
    starts_[base] = start;
    start += bwt_[base].Rank(text_size_);
  }
  key_starts_.Build();
}

EndFmIndex::Interval EndFmIndex::Extend(Interval range, unsigned base) const {
  if (range.begin == range.end) return {0, 0};
  return {starts_[base] + bwt_[base].Rank(range.begin),
          starts_[base] + bwt_[base].Rank(range.end)};
}

void EndFmIndex::BuildSeeds(int threads) {
  seed_length_ = std::min(length_, 8U);
  seeds_.push_back({0, text_size_});
  for (unsigned length = 1; length <= seed_length_; ++length) {
    size_t suffixes = seeds_.size();
    std::vector<Interval> next(suffixes * 4);
#pragma omp parallel for schedule(static) num_threads(threads) if(next.size() >= 256)
    for (size_t code = 0; code < next.size(); ++code)
      next[code] = Extend(seeds_[code % suffixes], code / suffixes);
    seeds_.swap(next);
  }
}

#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 6 && defined(__x86_64__)
__attribute__((target_clones("popcnt", "default")))
#endif
size_t EndFmIndex::Find(const uint64_t *words) const {
  if (keys_ == 0) return Missing();
  size_t code = 0;
  for (unsigned i = length_ - seed_length_; i < length_; ++i)
    code = (code << 2) | ((words[i / 32] >> (62 - 2 * (i % 32))) & 3);
  Interval range = seeds_[code];
  for (unsigned i = length_ - seed_length_; i > 0 && range.begin < range.end; --i) {
    unsigned base = (words[(i - 1) / 32] >> (62 - 2 * ((i - 1) % 32))) & 3;
    range = Extend(range, base);
  }
  return range.begin < range.end ? key_starts_.Rank(range.begin) : Missing();
}

size_t EndFmIndex::bytes() const {
  size_t total = text_.capacity() + seeds_.capacity() * sizeof(Interval) + key_starts_.bytes();
  for (const auto &rank : bwt_) total += rank.bytes();
  return total;
}
