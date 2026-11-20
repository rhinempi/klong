#include "end_anchors.h"
#include "platform.h"

#include <algorithm>
#include <exception>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_map>
#include <utility>
#include "edlib.h"
#include "assembly/anchor_paths.h"
#include "definitions.h"

namespace native {
namespace {
std::string field(const std::string &header, const std::string &name) {
  std::istringstream fields(header);
  std::string token;
  fields >> token;
  while (fields >> token)
    if (token.compare(0, name.size(), name) == 0) return token.substr(name.size());
  return "";
}

std::string reverse_complement(std::string seq) {
  std::reverse(seq.begin(), seq.end());
  for (auto &base : seq) {
    switch (base) {
      case 'A': base = 'T'; break;
      case 'C': base = 'G'; break;
      case 'G': base = 'C'; break;
      case 'T': base = 'A'; break;
    }
  }
  return seq;
}

std::string anchor(const std::string &header, const std::string &seq, size_t length) {
  const auto flag = field(header, "flag=");
  if (!flag.empty() && (std::stoul(flag) & contig_flag::kLoop)) return "";
  return seq.size() <= 2 * length ? seq : seq.substr(0, length) + seq.substr(seq.size() - length);
}

class FastaOutput {
 public:
  explicit FastaOutput(const std::string &path) : path_(path), file_(path + ".tmp") {
    if (!file_) throw Error("Cannot open FASTA output: " + path);
  }
  void write(const std::string &header, const std::string &seq) {
    file_ << header << '\n' << seq << '\n';
    ++records_;
    bases_ += seq.size();
  }
  void finish() {
    file_.close();
    if (!file_) throw Error("Cannot write FASTA output: " + path_);
    replace(path_ + ".tmp", path_);
    atomic_text(path_ + ".info", std::to_string(records_) + " " + std::to_string(bases_) + "\n");
  }
 private:
  std::string path_;
  std::ofstream file_;
  size_t records_{0}, bases_{0};
};

std::string graph_header(size_t id, const std::string &source_header) {
  const auto mul = field(source_header, "multi=");
  if (mul.empty()) throw Error("Missing contig multiplicity: " + source_header);
  return ">anchor_" + std::to_string(id) + " flag=0 multi=" + mul;
}

constexpr size_t kSeedLength = 16;
constexpr size_t kShortSeedLength = 12;
constexpr size_t kMinFragmentLength = 24;
constexpr size_t kBasesPerWord = 32;
unsigned base_code(char base) {
  switch (base) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 4;
  }
}

class PackedDna {
 public:
  explicit PackedDna(const std::string &seq)
      : length(seq.size()), words_((length + kBasesPerWord - 1) / kBasesPerWord, 0) {
    for (size_t i = 0; i < length; ++i) {
      auto code = base_code(seq[i]);
      if (code == 4) invalid.push_back(i);
      else words_[i / kBasesPerWord] |= uint64_t{code} << (2 * (i % kBasesPerWord));
    }
  }

  // Extract up to 32 bases at an arbitrary offset, including word boundaries.
  uint64_t word(size_t start, size_t count) const {
    const size_t index = start / kBasesPerWord, shift = 2 * (start % kBasesPerWord);
    uint64_t value = words_[index] >> shift;
    if (shift + 2 * count > 64) value |= words_[index + 1] << (64 - shift);
    if (count != kBasesPerWord) value &= (uint64_t{1} << (2 * count)) - 1;
    return value;
  }

  bool matches(size_t start, const PackedDna &query) const {
    for (size_t i = 0; i < query.length; i += kBasesPerWord)
      if (word(start + i, std::min(kBasesPerWord, query.length - i)) != query.words_[i / kBasesPerWord])
        return false;
    return true;
  }

  unsigned mismatches(size_t start, const PackedDna &query, size_t offset, size_t len) const {
    unsigned count = 0;
    for (size_t i = 0; i < len; i += kBasesPerWord) {
      size_t n = std::min(kBasesPerWord, len - i);
      uint64_t diff = word(start + i, n) ^ query.word(offset + i, n);
      count += __builtin_popcountll((diff | (diff >> 1)) & 0x5555555555555555ULL);
    }
    return count;
  }

  size_t length;
  std::vector<size_t> invalid;
 private:
  std::vector<uint64_t> words_;
};

constexpr size_t kMissing = std::numeric_limits<size_t>::max();
struct Placement {
  unsigned matches{0};
  unsigned mismatches{std::numeric_limits<unsigned>::max()};
  size_t preferred_unitig{kMissing}, preferred_hit{kMissing};
  size_t chain{kMissing}, member{0};
};
struct Query {
  size_t original;
  bool reverse;
  PackedDna sequence;
  bool stitched;
};
enum class CandidateKind { All, Threaded, RawShort };

bool eligible(const Query &query, CandidateKind kind) {
  return kind == CandidateKind::All || query.stitched == (kind == CandidateKind::Threaded);
}
struct Hit {
  size_t start, query;
  unsigned mismatches;
  unsigned span, columns, matches;
  Hit(size_t start_ = 0, size_t query_ = kMissing, unsigned mismatches_ = 0,
      unsigned span_ = 0, unsigned columns_ = 0, unsigned matches_ = 0)
      : start(start_), query(query_), mismatches(mismatches_), span(span_), columns(columns_), matches(matches_) {}
};
struct ShortSeed {
  size_t query, offset;
};
struct Alignment {
  size_t query, anchor_start, unitig_start, length;
  unsigned mismatches;
  bool exact_ends;
  unsigned span{0}, columns{0}, matches{0};
  Alignment(size_t q = 0, size_t a = 0, size_t u = 0, size_t len = 0, unsigned mm = 0, bool ends = false)
      : query(q), anchor_start(a), unitig_start(u), length(len), mismatches(mm), exact_ends(ends) {}
};
struct ReportHit {
  Alignment alignment;
  std::string unitig;
};
struct AnchorReport {
  ReportHit full, fragment;
};
struct MatchedUnitig {
  std::string sequence;
  std::vector<Hit> hits;
  CandidateKind kind;
  size_t owner, owner_start;
  std::string name;
};
struct HitLocation { size_t unitig, hit; };

size_t target_span(const Hit &hit, const std::vector<Query> &queries) {
  return hit.span ? hit.span : queries[hit.query].sequence.length;
}
struct Member {
  size_t original;
  bool reverse;
  std::string sequence;
};
struct Chain {
  size_t representative{kMissing};
  std::vector<Member> members;
  std::vector<std::string> bridges;
  std::string left, right;
};
using SeedIndex = std::unordered_map<uint32_t, std::vector<size_t>>;
using ShortSeedIndex = std::unordered_map<uint32_t, std::vector<ShortSeed>>;

std::vector<Hit> scan_unitig(const std::string &seq, const SeedIndex &seeds,
                             const std::vector<Query> &queries, CandidateKind kind) {
  const PackedDna packed(seq);
  std::vector<Hit> hits;
  auto invalid = packed.invalid.begin();
  for (size_t i = 0; i + kSeedLength <= seq.size(); ++i) {
    while (invalid != packed.invalid.end() && *invalid < i) ++invalid;
    const size_t valid_until = invalid == packed.invalid.end() ? seq.size() : *invalid;
    if (valid_until - i < kSeedLength) continue;
    auto it = seeds.find(packed.word(i, kSeedLength));
    if (it == seeds.end()) continue;
    for (size_t id : it->second) {
      const auto &query = queries[id];
      if (kind == CandidateKind::Threaded && !query.stitched) continue;
      if (query.sequence.length <= valid_until - i && packed.matches(i, query.sequence))
        hits.push_back({i, id, 0});
    }
  }
  return hits;
}

std::vector<Alignment> scan_short_seeds(const std::string &seq, const ShortSeedIndex &seeds,
                                       const std::vector<Query> &queries, size_t k, CandidateKind kind) {
  const PackedDna packed(seq);
  std::vector<std::pair<std::ptrdiff_t, size_t>> offsets;
  auto invalid = packed.invalid.begin();
  for (size_t pos = 0; pos + kShortSeedLength <= seq.size(); ++pos) {
    while (invalid != packed.invalid.end() && *invalid < pos) ++invalid;
    if (invalid != packed.invalid.end() && *invalid - pos < kShortSeedLength) continue;
    auto it = seeds.find(packed.word(pos, kShortSeedLength));
    if (it == seeds.end()) continue;
    for (const auto &seed : it->second) {
      if (kind == CandidateKind::Threaded && !queries[seed.query].stitched) continue;
      offsets.emplace_back(static_cast<std::ptrdiff_t>(pos) - static_cast<std::ptrdiff_t>(seed.offset), seed.query);
    }
  }
  std::sort(offsets.begin(), offsets.end());
  offsets.erase(std::unique(offsets.begin(), offsets.end()), offsets.end());

  std::vector<Alignment> hits;
  for (const auto &candidate : offsets) {
    const auto &anchor = queries[candidate.second].sequence;
    size_t start = candidate.first < 0 ? 0 : candidate.first;
    size_t offset = candidate.first < 0 ? -candidate.first : 0;
    size_t len = std::min(anchor.length - offset, seq.size() - start);
    if (len < std::min(kMinFragmentLength, anchor.length)) continue;
    if (!packed.invalid.empty()) {
      auto next_invalid = std::lower_bound(packed.invalid.begin(), packed.invalid.end(), start);
      if (next_invalid != packed.invalid.end() && *next_invalid < start + len) continue;
    }
    unsigned mm = packed.mismatches(start, anchor, offset, len);
    if (mm > len / 10) continue;
    const bool exact_ends = len == anchor.length &&
        (mm == 0 || (packed.mismatches(start, anchor, 0, k) == 0 &&
         packed.mismatches(start + len - k, anchor, len - k, k) == 0));
    hits.push_back({candidate.second, offset, start, len, mm, exact_ends});
  }
  return hits;
}

bool better_alignment(const Alignment &candidate, const Alignment &best, bool full) {
  if (!best.length) return true;
  if (full) return candidate.mismatches < best.mismatches;
  // Favor a long, high-identity fragment rather than a perfect but very short seed.
  size_t score = candidate.length - 2 * candidate.mismatches;
  size_t best_score = best.length - 2 * best.mismatches;
  return score > best_score || (score == best_score && candidate.length > best.length);
}

bool contains_export(const MatchedUnitig &outer, const Hit &outer_hit,
                     const MatchedUnitig &inner, const Hit &inner_hit,
                     const std::vector<Query> &queries) {
  const auto &query = queries[outer_hit.query];
  const auto &other = queries[inner_hit.query];
  const size_t len = query.sequence.length;
  auto left = [&](const MatchedUnitig &u, const Hit &h, bool reverse) {
    return reverse ? u.sequence.size() - h.start - len : h.start;
  };
  if (outer.sequence.size() < inner.sequence.size()) return false;
  size_t outer_left = left(outer, outer_hit, query.reverse);
  size_t inner_left = left(inner, inner_hit, other.reverse);
  if (outer_left < inner_left) return false;
  size_t offset = outer_left - inner_left;
  if (offset > outer.sequence.size() - inner.sequence.size()) return false;
  if (query.reverse) offset = outer.sequence.size() - inner.sequence.size() - offset;
  if (query.reverse == other.reverse)
    return outer.sequence.compare(offset, inner.sequence.size(), inner.sequence) == 0;
  for (size_t i = 0; i < inner.sequence.size(); ++i) {
    if (base_code(outer.sequence[offset + i]) + base_code(inner.sequence[inner.sequence.size() - 1 - i]) != 3)
      return false;
  }
  return true;
}

void coalesce_threaded_placements(const std::vector<MatchedUnitig> &matched,
                                  const std::vector<Query> &queries,
                                  std::vector<Placement> &placements) {
  struct Repeated {
    HitLocation owner{kMissing, 0}, covering{kMissing, 0};
    size_t last_unitig{kMissing};
    bool conflict{false};
  };
  std::vector<Repeated> repeated(placements.size());
  for (size_t u = 0; u < matched.size(); ++u) {
    const auto &unitig = matched[u];
    if (unitig.kind != CandidateKind::Threaded || unitig.owner == kMissing ||
        placements[unitig.owner].matches < 2) continue;
    auto &r = repeated[unitig.owner];
    for (size_t h = 0; h < unitig.hits.size(); ++h) {
      const auto &hit = unitig.hits[h];
      const auto &query = queries[hit.query];
      if (hit.mismatches || query.original != unitig.owner || query.reverse || hit.start != unitig.owner_start)
        continue;
      if (r.owner.unitig != kMissing) r.conflict = true;
      r.owner = r.covering = {u, h};
    }
  }
  for (size_t u = 0; u < matched.size(); ++u) {
    const auto &unitig = matched[u];
    if (unitig.kind != CandidateKind::Threaded) continue;
    for (size_t h = 0; h < unitig.hits.size(); ++h) {
      const auto &hit = unitig.hits[h];
      if (hit.mismatches) continue;
      auto &r = repeated[queries[hit.query].original];
      if (r.owner.unitig == kMissing || r.conflict) continue;
      if (r.last_unitig == u) {
        r.conflict = true;
        continue;
      }
      r.last_unitig = u;
      if (unitig.sequence.size() > matched[r.covering.unitig].sequence.size()) r.covering = {u, h};
    }
  }
  for (const auto &unitig : matched) {
    if (unitig.kind != CandidateKind::Threaded) continue;
    for (const auto &hit : unitig.hits) {
      if (hit.mismatches) continue;
      auto &r = repeated[queries[hit.query].original];
      if (r.owner.unitig == kMissing || r.conflict) continue;
      const auto &covering = matched[r.covering.unitig];
      const auto &reference = covering.hits[r.covering.hit];
      // The longest copy must cover every other one, regardless of input order.
      if (!contains_export(covering, reference, unitig, hit, queries)) r.conflict = true;
    }
  }
  for (size_t i = 0; i < repeated.size(); ++i) {
    if (repeated[i].owner.unitig == kMissing || repeated[i].conflict) continue;
    placements[i].matches = 1;
    placements[i].preferred_unitig = repeated[i].owner.unitig;
  }
}

bool prefer_target(const MatchedUnitig &target, const Hit &hit, const MatchedUnitig &old,
                   const Hit &old_hit, const std::vector<Query> &queries) {
  if (target.sequence.size() != old.sequence.size()) return target.sequence.size() > old.sequence.size();
  int cmp = target.sequence.compare(old.sequence);
  if (cmp != 0) return cmp < 0;
  if (hit.start != old_hit.start) return hit.start < old_hit.start;
  return queries[hit.query].reverse < queries[old_hit.query].reverse;
}

void find_gapped_hits(std::vector<MatchedUnitig> &targets, const std::vector<Query> &queries,
                      const std::vector<std::string> &anchor_sequences, std::vector<Placement> &placements,
                      int min_identity, int threads) {
  struct BestAlignment {
    size_t target{kMissing};
    Hit hit;
    std::exception_ptr error;
  };
  struct Result {
    EdlibAlignResult value;
    ~Result() { edlibFreeAlignResult(value); }
  };
  for (const auto &target : targets) {
    if (target.sequence.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
      throw Error("Unitig too long for Edlib anchor alignment");
  }
  std::vector<BestAlignment> results(queries.size() / 2);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
  for (size_t i = 0; i < results.size(); ++i) {
    const size_t original = queries[2 * i].original;
    if (placements[original].matches) continue;  // An exact hit already scores 100%.
    auto &best = results[i];
    try {
      for (size_t q : {2 * i, 2 * i + 1}) {
        const auto &anchor = anchor_sequences[q];
        unsigned floor_matches = min_identity, floor_columns = 100;
        if (best.target != kMissing) {
          floor_matches = best.hit.matches;
          floor_columns = best.hit.columns;
        }
        for (size_t t = 0; t < targets.size(); ++t) {
          const auto &target = targets[t];
          if (target.sequence.size() < (anchor.size() * floor_matches + floor_columns - 1) / floor_columns)
            continue;
          // With d edits, identity cannot exceed query length / (query length + d).
          const int edits = std::min(anchor.size(), anchor.size() * (floor_columns - floor_matches) / floor_matches);
          Result alignment{edlibAlign(anchor.data(), anchor.size(), target.sequence.data(), target.sequence.size(),
                                      edlibNewAlignConfig(edits, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0))};
          const auto &aln = alignment.value;
          if (aln.status != EDLIB_STATUS_OK) throw Error("Edlib failed aligning a flanking anchor");
          if (aln.editDistance < 0 || !aln.alignmentLength || aln.endLocations[0] < aln.startLocations[0]) continue;
          unsigned matches = std::count(aln.alignment, aln.alignment + aln.alignmentLength, EDLIB_EDOP_MATCH);
          unsigned substitutions = std::count(aln.alignment, aln.alignment + aln.alignmentLength, EDLIB_EDOP_MISMATCH);
          Hit hit(aln.startLocations[0], q, substitutions, aln.endLocations[0] - aln.startLocations[0] + 1,
                  aln.alignmentLength, matches);
          auto score = uint64_t{matches} * floor_columns;
          auto floor = uint64_t{floor_matches} * hit.columns;
          if (score < floor || (score == floor && best.target != kMissing &&
              !prefer_target(target, hit, targets[best.target], best.hit, queries))) continue;
          best.target = t;
          best.hit = hit;
          floor_matches = matches;
          floor_columns = hit.columns;
        }
      }
    } catch (...) {
      best.error = std::current_exception();
    }
  }
  for (const auto &result : results) {
    if (result.error) std::rethrow_exception(result.error);
    if (result.target == kMissing) continue;
    placements[queries[result.hit.query].original].mismatches = result.hit.mismatches;
    targets[result.target].hits.push_back(result.hit);
  }
  for (auto &target : targets) std::sort(target.hits.begin(), target.hits.end(), [](const Hit &a, const Hit &b) {
    return a.start != b.start ? a.start < b.start : a.query < b.query;
  });
}

void choose_best_hits(const std::vector<MatchedUnitig> &matched, const std::vector<Query> &queries,
                      std::vector<Placement> &placements, std::vector<AnchorReport> &reports) {
  for (size_t u = 0; u < matched.size(); ++u) {
    const auto &unitig = matched[u];
    for (size_t h = 0; h < unitig.hits.size(); ++h) {
      const auto &hit = unitig.hits[h];
      auto &placement = placements[queries[hit.query].original];
      if (hit.mismatches != placement.mismatches) continue;
      if (placement.preferred_unitig != kMissing) {
        const auto &old = matched[placement.preferred_unitig];
        const auto &old_hit = old.hits[placement.preferred_hit];
        if (!prefer_target(unitig, hit, old, old_hit, queries)) continue;
      }
      placement.matches = 1;
      placement.preferred_unitig = u;
      placement.preferred_hit = h;
      if (!reports.empty()) {
        reports[queries[hit.query].original].full = {
            {hit.query, 0, hit.start, queries[hit.query].sequence.length, hit.mismatches, true}, unitig.name};
        auto &alignment = reports[queries[hit.query].original].full.alignment;
        alignment.span = hit.span;
        alignment.columns = hit.columns;
        alignment.matches = hit.matches;
      }
    }
  }
}

void select_best_sites(MatchedUnitig &unitig, size_t unitig_id, const std::vector<Query> &queries,
                       const std::vector<Placement> &placements) {
  const auto &hits = unitig.hits;
  auto len = [&](size_t i) { return target_span(hits[i], queries); };
  auto end = [&](size_t i) { return hits[i].start + len(i); };
  std::vector<size_t> ranked;
  for (size_t i = 0; i < hits.size(); ++i) {
    const auto &p = placements[queries[hits[i].query].original];
    if (p.preferred_unitig == unitig_id && p.preferred_hit == i) ranked.push_back(i);
  }
  std::sort(ranked.begin(), ranked.end(), [&](size_t a, size_t b) {
    if (hits[a].columns || hits[b].columns) {
      auto cols = [](const Hit &hit) { return hit.columns ? hit.columns : 1u; };
      auto matches = [](const Hit &hit) { return hit.columns ? hit.matches : 1u; };
      auto lhs = uint64_t{matches(hits[a])} * cols(hits[b]);
      auto rhs = uint64_t{matches(hits[b])} * cols(hits[a]);
      if (lhs != rhs) return lhs > rhs;
    }
    if (hits[a].mismatches != hits[b].mismatches) return hits[a].mismatches < hits[b].mismatches;
    if (len(a) != len(b)) return len(a) > len(b);
    return a < b;
  });

  std::set<size_t> chosen;
  for (size_t i : ranked) {
    auto next = chosen.lower_bound(i);
    if ((next == chosen.end() || end(i) <= hits[*next].start) &&
        (next == chosen.begin() || end(*std::prev(next)) <= hits[i].start)) chosen.insert(next, i);
  }
  std::vector<Hit> sites;
  for (size_t i : chosen) sites.push_back(hits[i]);
  for (size_t i = 0; !chosen.empty() && i < hits.size(); ++i) {
    // Ignore alternative placements as hits, but don't copy through a known
    // artificial join unless the copied span is replaced with an original.
    if (!queries[hits[i].query].stitched || hits[i].columns || hits[i].mismatches || chosen.count(i)) continue;
    size_t join = hits[i].start + len(i) / 2;
    auto next = std::lower_bound(sites.begin(), sites.begin() + chosen.size(), join,
                                 [](const Hit &hit, size_t p) { return hit.start < p; });
    if (next != sites.begin()) {
      const auto &prev = *std::prev(next);
      if (join < prev.start + target_span(prev, queries)) continue;
    }
    sites.push_back({join, kMissing, 0});
  }
  std::sort(sites.begin(), sites.end(), [](const Hit &a, const Hit &b) {
    if (a.start != b.start) return a.start < b.start;
    return a.query == kMissing && b.query != kMissing;
  });
  unitig.hits = std::move(sites);
}

void plan_chains(const MatchedUnitig &unitig, size_t unitig_id, const std::vector<Query> &queries,
                 std::vector<Placement> &placements, std::vector<Chain> &chains) {
  const auto &seq = unitig.sequence;
  const auto &hits = unitig.hits;
  auto end = [&](size_t i) {
    return hits[i].start + (hits[i].query == kMissing ? 0 : target_span(hits[i], queries));
  };
  auto next_start = [&](size_t i) { return i + 1 < hits.size() ? hits[i + 1].start : seq.size(); };
  std::vector<size_t> earlier_end(hits.size(), 0);
  for (size_t i = 1; i < hits.size(); ++i)
    earlier_end[i] = std::max(earlier_end[i - 1], end(i - 1));
  auto usable = [&](size_t i) {
    if (hits[i].query == kMissing) return false;
    const auto &query = queries[hits[i].query];
    const auto &placement = placements[query.original];
    return eligible(query, unitig.kind) && placement.matches == 1 &&
           (placement.preferred_unitig == kMissing || placement.preferred_unitig == unitig_id) &&
           hits[i].start >= earlier_end[i] && end(i) <= next_start(i);
  };

  for (size_t i = 0; i < hits.size(); ++i) {
    if (!usable(i)) continue;
    const size_t first = i;
    Chain chain;
    // Unusable matches are barriers, not anchors that can be skipped over.
    for (;;) {
      const auto &query = queries[hits[i].query];
      auto &placement = placements[query.original];
      placement.chain = chains.size();
      placement.member = chain.members.size();
      chain.members.push_back({query.original, query.reverse, {}});
      chain.representative = std::min(chain.representative, query.original);
      if (i + 1 == hits.size() || !usable(i + 1) ||
          !std::all_of(seq.begin() + end(i), seq.begin() + next_start(i),
                       [](char c) { return base_code(c) < 4; })) break;
      chain.bridges.push_back(seq.substr(end(i), next_start(i) - end(i)));
      ++i;
    }
    size_t left = hits[first].start, right = end(i);
    while (left > earlier_end[first] && base_code(seq[left - 1]) < 4) --left;
    while (right < next_start(i) && base_code(seq[right]) < 4) ++right;
    chain.left = seq.substr(left, hits[first].start - left);
    chain.right = seq.substr(end(i), right - end(i));
    if (chain.members.size() == 1 && chain.members.front().reverse) {
      std::swap(chain.left, chain.right);
      chain.left = reverse_complement(std::move(chain.left));
      chain.right = reverse_complement(std::move(chain.right));
    }
    chains.push_back(std::move(chain));
  }
}

}  // namespace

void write_end_anchors(const std::string &source, const std::string &output, size_t end_length) {
  FastaOutput out(output);
  size_t id = 0;
  read_fasta(source, [&](const std::string &header, const std::string &seq) {
    const auto truncated = anchor(header, seq, end_length);
    if (!truncated.empty()) out.write(graph_header(id++, header), truncated);
  });
  out.finish();
}

void write_outward_consensus(const std::string &source, const std::string &output, size_t end_length) {
  FastaOutput out(output);
  size_t id = 0;
  read_fasta(source, [&](const std::string &header, const std::string &seq) {
    size_t left = std::stoul(field(header, "left_extension="));
    size_t right = std::stoul(field(header, "right_extension="));
    if (left > seq.size() || right > seq.size() - left) throw Error("Invalid consensus lengths: " + header);
    size_t context = std::min(end_length, seq.size() - left - right);
    if (left) out.write(graph_header(id++, header), seq.substr(0, left + context));
    if (right) out.write(graph_header(id++, header), seq.substr(seq.size() - right - context));
  });
  out.finish();
}

AnchorRestoreStats restore_end_anchors(const std::string &source, const std::string &unitigs,
                                       const std::string &output, size_t end_length, size_t k, int threads,
                                       int max_mismatches, const std::string &report,
                                       const std::string &short_unitigs, bool best_hit, int min_identity) {
  if (!k || k > kMaxK || k > end_length || end_length > kFlankingAnchorEndLength || threads < 1 ||
      max_mismatches < -1 || max_mismatches > 8)
    throw Error("Invalid end anchor length, graph k, thread count, or mismatch limit");
  if (best_hit && !short_unitigs.empty()) throw Error("Best-hit restoration requires raw unitigs only");
  const bool gapped = min_identity > 0;
  if (min_identity < 0 || min_identity > 100 || (gapped && (!best_hit || max_mismatches >= 0)))
    throw Error("Gapped anchor alignment requires raw best-hit restoration without the short-seed option");
  const bool approximate = max_mismatches > 0;
  std::vector<Placement> placements;
  std::vector<AnchorReport> reports;
  std::vector<Query> queries;
  std::vector<std::string> anchor_sequences;
  SeedIndex seeds;
  ShortSeedIndex short_seeds;
  std::unordered_map<std::string, size_t> threaded_owners;
  read_fasta(source, [&](const std::string &header, const std::string &seq) {
    size_t original = placements.size();
    placements.emplace_back();
    if (!report.empty()) reports.emplace_back();
    auto key = anchor(header, seq, end_length);
    if (key.size() <= k || key.size() < kSeedLength) return;
    if (!short_unitigs.empty() && seq.size() > 2 * end_length) {
      auto entry = threaded_owners.emplace(header.substr(1, header.find_first_of(" \t") - 1), original);
      if (!entry.second) entry.first->second = kMissing;
    }
    for (bool reverse : {false, true}) {
      PackedDna packed(key);
      if (!packed.invalid.empty()) return;
      if (!approximate) seeds[packed.word(0, kSeedLength)].push_back(queries.size());
      if (!gapped && (approximate || !report.empty())) {
        // Overlapping tiles also find partial anchors that start between two tiles.
        for (size_t i = 0; i + kShortSeedLength <= packed.length; i += kShortSeedLength / 2)
          short_seeds[packed.word(i, kShortSeedLength)].push_back({queries.size(), i});
      }
      queries.push_back({original, reverse, std::move(packed), seq.size() > 2 * end_length});
      if (gapped) anchor_sequences.push_back(key);
      if (!reverse) key = reverse_complement(std::move(key));
    }
  });

  std::vector<std::string> batch;
  std::vector<std::string> unitig_names;
  std::vector<std::pair<size_t, size_t>> batch_owners;
  std::vector<MatchedUnitig> matched;
  size_t batch_bases = 0;
  auto flush = [&](CandidateKind kind) {
    if (batch.empty()) return;
    std::vector<std::vector<Hit>> results(batch.size());
    std::vector<std::vector<Alignment>> alignments;
    if (!short_seeds.empty()) alignments.resize(batch.size());
    std::vector<std::exception_ptr> errors(batch.size());
    const int workers = std::min(threads, static_cast<int>(batch.size()));
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers) if(workers > 1)
    for (size_t i = 0; i < batch.size(); ++i) {
      try {
        if (!seeds.empty()) results[i] = scan_unitig(batch[i], seeds, queries, kind);
        if (!alignments.empty()) alignments[i] = scan_short_seeds(batch[i], short_seeds, queries, k, kind);
      } catch (...) {
        errors[i] = std::current_exception();
      }
    }
    for (const auto &error : errors) if (error) std::rethrow_exception(error);
    // Track the best mismatch count and cap the number of equally good hits.
    for (size_t i = 0; i < results.size(); ++i) {
      if (!alignments.empty()) {
        for (const auto &candidate : alignments[i]) {
          const auto &query = queries[candidate.query];
          const bool full = candidate.length == query.sequence.length;
          if (!report.empty() && eligible(query, kind)) {
            auto &summary = reports[query.original];
            auto &best = full ? summary.full : summary.fragment;
            if (better_alignment(candidate, best.alignment, full))
              best = {candidate, unitig_names[i]};
          }
          // Outward flanks must attach to an unchanged terminal k-mer on either side.
          if (approximate && candidate.exact_ends && candidate.mismatches <= static_cast<unsigned>(max_mismatches))
            results[i].push_back({candidate.unitig_start, candidate.query, candidate.mismatches});
        }
      }
      bool usable_candidate = false;
      for (const auto &hit : results[i]) {
        const auto &query = queries[hit.query];
        // In raw short candidates, a matched stitched anchor blocks flank copying.
        if (!eligible(query, kind)) continue;
        usable_candidate = true;
        auto &placement = placements[query.original];
        if (hit.mismatches < placement.mismatches) {
          placement.mismatches = hit.mismatches;
          placement.matches = 1;
        } else if (hit.mismatches == placement.mismatches && placement.matches < 2) ++placement.matches;
      }
      if (usable_candidate || gapped) {
        auto owner = kind == CandidateKind::Threaded ? batch_owners[i] : std::make_pair(kMissing, kMissing);
        matched.push_back({std::move(batch[i]), std::move(results[i]), kind, owner.first, owner.second,
                           best_hit && !report.empty() ? std::move(unitig_names[i]) : std::string()});
      }
    }
    batch.clear();
    unitig_names.clear();
    batch_owners.clear();
    batch_bases = 0;
  };
  const auto scan_candidates = [&](const std::string &path, CandidateKind kind) {
    read_fasta(path, [&](const std::string &header, const std::string &seq) {
      if ((seeds.empty() && short_seeds.empty()) || seq.empty() || (!gapped && seq.size() < kShortSeedLength)) return;
      if (batch.size() >= 256 || batch_bases + seq.size() > (4ULL << 20)) flush(kind);
      batch.emplace_back(seq);
      if (!report.empty() || kind == CandidateKind::Threaded) {
        auto name = header.substr(1, header.find_first_of(" \t") - 1);
        if (kind == CandidateKind::Threaded) {
          auto owner = threaded_owners.find(name);
          auto start = field(header, "left_extension=");
          batch_owners.emplace_back(owner == threaded_owners.end() || start.empty() ? kMissing : owner->second,
                                    start.empty() ? kMissing : std::stoul(start));
        }
        if (!report.empty()) unitig_names.push_back(std::move(name));
      }
      batch_bases += seq.size();
    });
    flush(kind);
  };
  scan_candidates(unitigs, short_unitigs.empty() ? CandidateKind::All : CandidateKind::Threaded);
  if (!short_unitigs.empty()) scan_candidates(short_unitigs, CandidateKind::RawShort);
  if (!short_unitigs.empty()) coalesce_threaded_placements(matched, queries, placements);
  if (gapped && min_identity < 100) find_gapped_hits(matched, queries, anchor_sequences, placements, min_identity, threads);
  if (best_hit) choose_best_hits(matched, queries, placements, reports);

  std::vector<Chain> chains;
  for (size_t u = 0; u < matched.size(); ++u) {
    auto &unitig = matched[u];
    if (best_hit) {
      select_best_sites(unitig, u, queries, placements);
    } else if (approximate) {
      unitig.hits.erase(std::remove_if(unitig.hits.begin(), unitig.hits.end(), [&](const Hit &hit) {
        const auto &query = queries[hit.query];
        return eligible(query, unitig.kind) && hit.mismatches != placements[query.original].mismatches;
      }), unitig.hits.end());
    }
    plan_chains(unitig, u, queries, placements, chains);
  }
  matched.clear();
  AnchorRestoreStats stats;
  for (const auto &chain : chains) stats.joined += chain.members.size() - 1;
  if (stats.joined) {
    // Only originals used in multi-anchor chains need to be held for reordering.
    size_t id = 0;
    read_fasta(source, [&](const std::string &, const std::string &seq) {
      const auto &placement = placements[id++];
      if (placement.chain == kMissing) return;
      auto &chain = chains[placement.chain];
      if (chain.members.size() == 1) return;
      auto &member = chain.members[placement.member];
      member.sequence = member.reverse ? reverse_complement(seq) : seq;
    });
  }
  FastaOutput out(output);
  std::ofstream details;
  if (!report.empty()) {
    details.open(report + ".tmp");
    if (!details) throw Error("Cannot open anchor match report: " + report);
    details << "source_contig\tanchor_nt\tmatch\tunitig\tstrand\tanchor_start\tunitig_start"
               "\taligned_nt\tmismatches\tunaligned_nt";
    if (gapped) details << "\ttarget_span_nt\talignment_columns\tmatched_nt\tgap_nt\tidentity_pct\trestored";
    details << '\n';
  }
  read_fasta(source, [&](const std::string &header, const std::string &seq) {
    size_t id = stats.contigs++;
    if (!report.empty()) {
      const auto &summary = reports[id];
      const auto name = header.substr(1, header.find_first_of(" \t") - 1);
      const auto length = anchor(header, seq, end_length).size();
      auto row = [&](const char *kind, const ReportHit &best) {
        const auto &hit = best.alignment;
        if (!hit.length) return;
        const bool reverse = queries[hit.query].reverse;
        const size_t begin = reverse ? length - hit.anchor_start - hit.length : hit.anchor_start;
        details << name << '\t' << length << '\t' << kind << '\t' << best.unitig << '\t'
                << (reverse ? '-' : '+') << '\t' << begin + 1 << '\t' << hit.unitig_start + 1
                << '\t' << hit.length << '\t' << hit.mismatches << '\t' << length - hit.length;
        if (gapped) {
          const auto columns = hit.columns ? hit.columns : hit.length;
          const auto matches = hit.columns ? hit.matches : hit.length - hit.mismatches;
          details << '\t' << (hit.span ? hit.span : hit.length) << '\t' << columns << '\t' << matches
                  << '\t' << columns - matches - hit.mismatches << '\t' << 100.0 * matches / columns
                  << '\t' << (placements[id].chain != kMissing);
        }
        details << '\n';
      };
      row("full", summary.full);
      row("fragment", summary.fragment);
      if (!summary.full.alignment.length && !summary.fragment.alignment.length) {
        details << name << '\t' << length << "\tno_hit\t.\t.\t.\t.\t0\t.\t" << length;
        if (gapped) details << "\t.\t.\t.\t.\t.\t0";
        details << '\n';
      }
    }
    const auto &placement = placements[id];
    if (placement.chain == kMissing) {
      out.write(header, seq);
      return;
    }
    const auto &chain = chains[placement.chain];
    if (id != chain.representative) return;
    size_t added = chain.left.size() + chain.right.size();
    for (const auto &bridge : chain.bridges) added += bridge.size();
    if (chain.members.size() == 1 && !added) {
      out.write(header, seq);
      return;
    }
    if (added) ++stats.extended;
    stats.bases_added += added;
    std::string restored;
    size_t length = added;
    if (chain.members.size() == 1) length += seq.size();
    else for (const auto &member : chain.members) length += member.sequence.size();
    restored.reserve(length);
    restored += chain.left;
    if (chain.members.size() == 1) restored += seq;
    else {
      for (size_t i = 0; i < chain.members.size(); ++i) {
        if (i) restored += chain.bridges[i - 1];
        restored += chain.members[i].sequence;
      }
    }
    restored += chain.right;
    std::string updated = header;
    auto len = updated.find(" len=");
    if (len != std::string::npos) {
      auto end = updated.find_first_of(" \t", len + 5);
      updated.replace(len + 5, end - (len + 5), std::to_string(restored.size()));
    } else updated += " len=" + std::to_string(restored.size());
    out.write(updated, restored);
  });
  out.finish();
  if (!report.empty()) {
    details.close();
    if (!details) throw Error("Cannot write anchor match report: " + report);
    replace(report + ".tmp", report);
  }
  return stats;
}
}  // namespace native
