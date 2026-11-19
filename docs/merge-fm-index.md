# Contig-End FM Index

Native reduction now retains only contigs longer than next-k and extracts their
edges directly, without short-contig matching or an overlap index. The native
`--merge-index` setting remains accepted and saved for compatibility, but does
not affect these long-only merges. The implementation and historical benchmarks
below apply to the standalone core's optional short-contig overlap mode.

The merge stage supports two interchangeable exact-overlap dictionaries:

```sh
megahit --merge-index fm [other assembly options]
megahit --merge-index hash [other assembly options]
```

Hash remains the default: measurements below do not justify a blanket switch to
FM for fixed-length exact overlap. The FM backend replaces the query hash tables
when selected, rather than adding an approximate or alignment-based lookup.
The setting is stored in native `options.json`. A resumed run uses its saved
options and existing checkpoints.

## FM Representation

Each query collection is encoded as DNA records separated by a non-DNA symbol,
with a unique final sentinel. The vendored libsais v2.12.1 library constructs a
suffix array using OpenMP. The implementation derives the BWT and four DNA rank
dictionaries, then discards the suffix array and original text. Separator symbols
prevent matches across contig boundaries. Equal ends share a single query ID.

An 8-base interval table accelerates backward search. Remaining query bases are
checked by exact FM interval updates. A bit-vector rank maps a matching interval
to its query ID without a hash table or retained full suffix array. The index is
immutable during the scan, so workers query it without locks. GCC/x86-64 builds
include safely dispatched POPCNT and generic lookup variants.

The FM index is exact. When `--merge-dedup bloom` is selected, the two existing
Bloom filters still make extension redundancy filtering probabilistic. Index
selection does not alter their hash policy or sizing. The independent
default [`--merge-dedup sort` backend](merge-sort-dedup.md) removes these filters entirely.

## Parallelism And Ordering

Suffix-array construction, BWT bit-vector construction, and seed-table generation
use the requested worker count where the workload is large enough. Read-edge
search uses the existing parallel shard ownership when there are enough files.
With fewer files than workers, base edges are first copied in original order;
then record-aligned ranges are scanned in parallel. This adds an input pass but
allows one large shard to use multiple search workers.

Matches carry resolved prefix/suffix IDs, avoiding repeated dictionary searches
during consolidation. Temporary match sequences use only active packed words
(8 bytes for a k29 edge, not 64), plus two IDs. Per-range storage is released as
soon as it is consolidated. Final query-hit vectors still contain fixed-width
keys and can grow large for highly repetitive ends; matching is not disk-backed.

Range results are consolidated in input order. Extension generation and emission
remain ordered and serial to preserve duplicate filtering and output order.
Diagnostic mode retains a serial scan for its insertion-order-dependent counts.
Bloom initialization also retains its existing merge behavior. This is not an
entirely parallel merge pipeline.

## Validation

The merge oracle suite exercises both backends through k255, duplicated ends,
palindromes, seed collisions, invalid bases, empty inputs, and several worker
counts. Binary output and diagnostic extension reports are compared byte-for-byte
with the previous executable when `MEGAHIT_TEST_BASELINE_CORE` is set. Native
workflow tests exercise both command-line backend selections and saved options.

## Benchmarks

The optional native `prefixbf --mode index-bench` mode measures dictionary build
and a parallel scan without allocating Bloom filters or writing edges. Scan time
includes input reading and packed prefix/suffix extraction, not just a lookup
instruction. Both backends must report identical key, orientation, and hit counts.

Three alternating runs per backend/thread count, using warm filesystem caches:

| Dataset | Threads | Hash Scan | FM Scan | Hash Build | FM Build |
| --- | ---: | ---: | ---: | ---: | ---: |
| E. coli, 4,579,364 edges, about 1,052 keys per end | 1 | 0.916 s | 0.727 s | 0.045 s | 0.051 s |
| Same E. coli inputs | 8 | 0.277 s | 0.434 s | 0.021 s | 0.028 s |
| Synthetic, 500,000 edges, 100,000 keys per end | 1 | 0.206 s | 0.315 s | 0.142 s | 0.314 s |
| Same synthetic inputs | 8 | 0.102 s | 0.137 s | 0.119 s | 0.282 s |

Values are medians. Lookup-only process peak RSS at eight threads was 14.7 MiB
(hash) versus 17.2 MiB (FM) for E. coli, and 24.1 versus 27.1 MiB for synthetic
inputs. Construction workspace contributes to RSS; these are not full merge
memory figures. A full E. coli FM merge still uses roughly 580 MiB, primarily
because the existing Bloom filters are unchanged. Small subsecond differences
are sensitive to CPU placement, cache, and system load; they are not predictions
for a marine metagenome.

FM improved the single-threaded E. coli scan but did not improve the measured
parallel scans or larger synthetic dictionary. Fixed-length integer hashing
remains competitive; BWT indexes are not inherently faster for this workload.

Reproduce the synthetic comparison with the optional Python test harness:

```sh
python3 tests/benchmark_merge_index.py --core build/megahit_core \
  --work-dir /tmp/merge-index-benchmark --queries 100000 --edges 500000 \
  --threads 1,8 --repeats 3
```

For real inputs, supply `--long-prefix PREFIX --contig SHORT_CONTIG_FASTA`
instead of generating synthetic fixtures. Results and per-run peak RSS are saved
in `results.json` under the work directory. Native production builds and runtime
still require no Python.
