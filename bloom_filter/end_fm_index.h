#ifndef MEGAHIT_END_FM_INDEX_H
#define MEGAHIT_END_FM_INDEX_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// Exact fixed-length DNA dictionary. Build once, then query concurrently.
class EndFmIndex {
 public:
  static size_t Missing() { return std::numeric_limits<size_t>::max(); }
  explicit EndFmIndex(unsigned length) : length_(length) {}
  void Add(const uint64_t *words);
  void Build(int threads);
  size_t Find(const uint64_t *words) const;
  size_t size() const { return keys_; }
  size_t bytes() const;

 private:
  struct RankBits {
    std::vector<uint64_t> bits, counts;
    void Init(size_t length);
    void Build();
    uint64_t Rank(size_t end) const;
    size_t bytes() const { return 8 * (bits.capacity() + counts.capacity()); }
  };
  struct Interval { uint64_t begin, end; };
  template <typename Index> void Finish(const std::vector<Index> &sa, int threads);
  Interval Extend(Interval range, unsigned base) const;
  void BuildSeeds(int threads);

  unsigned length_, seed_length_ = 0;
  size_t keys_ = 0, text_size_ = 0;
  std::vector<uint8_t> text_;
  std::array<RankBits, 4> bwt_;
  std::array<uint64_t, 4> starts_{};
  RankBits key_starts_;
  std::vector<Interval> seeds_;
};

#endif
