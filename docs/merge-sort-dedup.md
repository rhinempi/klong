# Exact Sort-Based Merge Deduplication

The current [native workflow](native-executable.md) uses exact sorting for
`reduced_edges` but no longer runs the edge merge described below.
`--merge-sort-memory` still controls reduced-edge sort buffers.

The previous forward-only workflow deferred intermediate merge deduplication until after
reduction. With the default 21,29,39,59 list, `k29/merged_edges` and
`k39/merged_edges` retain duplicates without Bloom filtering or sorting. Only
surviving records enter an exact external sorter and are written directly to
`kN/reduced_edges` before graph construction. Reverse complements share one
canonical packed key. Intermediate merges write read-counted and generated edges
to separate shards. When both sources supply the same surviving edge, reduction
deduplication keeps the original read count, even if it is below 2. Within a
source, duplicate survivors keep the maximum multiplicity, not the sum.

In that workflow, the last reduction k (normally k59) had no subsequent
reduction. Its `merged_edges` used the configured merge backend, exact sort
by default, with the original base-record preservation rule described below.
To select the standalone core's exact merge backend:

```sh
megahit --core prefixbf --mode merge --merge_dedup sort --sort_memory_mb 128 [edge and contig inputs]
```

`--merge-sort-memory` is the total record-buffer budget in MiB, shared by sorting
workers, not a per-worker budget. Its default is 128 MiB and its minimum is
0.0625 MiB. The budget is capped at one quarter of the smaller visible available
memory and workflow `-m` budget. Indexes, overlap matches, I/O buffers, merge
cursors, and runtime overhead are additional; this is not an enforced process
RSS limit. The standalone core's legacy short-contig overlap mode can still use
`--merge_index hash|fm`.

Bloom remains selectable in standalone merge mode with `--merge_dedup bloom`.
This option does not change the separate k-mer reduction Bloom filter.
The installed native executable still requires no Python, GNU sort, or other
external sorting program.

## Parallelism

Long-only contig extraction uses OpenMP workers and FASTA batches bounded to
8 MiB of sequence plus the final record. Packed candidates enter thread-private
sort buffers. Reduction inserts surviving records directly into per-worker
sort buffers, avoiding a separate unfiltered-library sorting pass.

Buffer sorting/flushing uses up to `-t` workers. Independent external merge
groups run concurrently. Final key merging uses disjoint key ranges, sampled
boundaries and binary seeks into each input run; equal keys never cross range
boundaries. Merge concurrency is capped at `min(-t, 8)` to bound open files and
I/O buffers (at most 32 input runs per merge worker). Reduced unique records are
emitted concurrently into separate edge shards; ordering and shard distribution
can differ between thread counts.

The terminal exact merge parallelizes key selection and the buffers for restoring
candidate generation order, but its final ordinal merge/output remains serial
to preserve historical output order. FASTA parsing and each individual
`std::sort` remain serial; separate worker buffers sort concurrently.

## Terminal Merge

1. Copy every read-counted base edge unchanged to its original output shard.
2. Collect its canonical packed key with ordinal zero into per-worker buffers.
3. Generate candidates by extracting packed (next-k+1)-mer windows from retained
   contigs longer than next-k. Store canonical keys with monotonically increasing
   generation ordinals. The standalone core also retains its optional
   short-contig exact-overlap candidate generator.
4. Sort full buffers into binary runs using `std::sort`, discarding equal-key
   duplicates within each run. Flush remaining worker buffers in parallel.
5. Merge runs by key and ordinal. Base ordinal zero takes precedence, suppressing
   any candidate already present in the input. Otherwise retain the earliest
   candidate for each key. Multi-pass merging opens at most 32 input runs per worker.
6. Sort accepted candidates by their retained ordinals and emit them in original
   generation order, with round-robin shard assignment. Standalone core merges
   assign generated edges multiplicity 50 unless overridden.

Records store only active 64-bit sequence words plus one ordinal: 16 bytes at
k29 and 24 bytes at k39/k59. Neither merge Bloom filter is allocated. Base input
multiplicities, orientations, and even pre-existing duplicate base records are
preserved; this backend makes newly added edges unique, rather than rewriting
the base library or aggregating its counts.

The external sorter is in `bloom_filter/exact_merge_dedup.cpp`. Temporary run
directories are siblings of the merged edge prefix, for example
`tmp/k29/reduced_edges/29.sort-XXXXXX` or `tmp/k59/merged_edges/59.sort-XXXXXX`.
They are removed on success and handled
exceptions. Abrupt termination can leave a workspace behind. Run read/write
failures are checked, and incomplete work does not publish new edge metadata.
Per-candidate provenance text reports are not supported with sort mode; the core
requires `--human_readable 0`, as the native workflow already supplies.

Sorting is exact, so Bloom false positives cannot discard candidates. Byte-for-
byte equality with Bloom is expected when its membership answers are correct,
but is not guaranteed on arbitrary data because Bloom can have false positives.

## Benchmarks

The paired-end E. coli smoke test used 8 threads, a 4,000,000,000-byte memory
budget, default phase lists, and `--flanking-consensus-only --keep-tmp-files`.
Its k29 merge retained 4,554,694 records without sorting. Reduction left 2,459
records for the sorter, which wrote 2,432 unique edges and removed 27 duplicates.
The sorter wrote 39,344 temporary bytes on a focused rerun. The k39 reduced
output contained 9,649 unique edges; the terminal k59 merge still used exact
deduplication. Full assembly completed in 60.1 seconds with 135 final contigs,
4,579,750 bp and N50 94,524 bp. This is an end-to-end smoke test, not a controlled
scaling or before/after benchmark.

The measurements below predate deferred reduction deduplication and parallel
key-range merging; they describe the earlier standalone merge implementation.

Complete merge processes were measured with three interleaved Bloom/sort runs,
the same hash overlap index and input edges, a 4,000,000,000-byte workflow budget,
and 128 MiB of sort record buffers. Times include setup, indexing, scanning,
candidate generation, exact sorting, output, and temporary cleanup. Peak memory
is process RSS, not filesystem page cache. Input was warm and temporary I/O was
buffered on the local filesystem; this is not a cold-disk or durable-write test.

Eight-worker medians:

| Merge | Bloom Time | Sort Time | Bloom RSS | Sort RSS | Sort Temp Writes |
| --- | ---: | ---: | ---: | ---: | ---: |
| E. coli k21 to k29 | 1.56 s | 1.01 s | 581.9 MiB | 87.7 MiB | 70.0 MiB |
| E. coli k29 to k39 | 1.51 s | 1.11 s | 583.0 MiB | 121.5 MiB | 105.4 MiB |
| E. coli k39 to k59 | 1.67 s | 1.29 s | 583.8 MiB | 123.2 MiB | 106.0 MiB |
| Synthetic, 500k base edges / 625k candidates | 1.27 s | 0.96 s | 122.9 MiB | 73.4 MiB | 26.7 MiB |

Temporary reads equaled the writes in these runs. These are total temporary bytes
written, not peak temporary disk occupancy; final edge output is additional.
Every output shard and metadata file had the same SHA-256 as the corresponding
Bloom run, including edge order and multiplicities.

One-worker medians, Bloom versus sort: k21-to-k29 4.21 versus 1.72 s; k29-to-k39
4.49 versus 2.03 s; k39-to-k59 5.52 versus 2.87 s; synthetic 1.76 versus 1.26 s.
The synthetic input produced 625,000 candidates and 623,740 new edges, exercising
a candidate-heavy case rather than only the relatively sparse E. coli extensions.

These results favor sorting over the current merge Bloom implementation, not
over every possible tuned Bloom design. Sorting adds disk traffic; inputs much
larger than RAM, smaller sort buffers, or slow/network scratch storage can require
more merge passes and change the result. Match-vector memory can still grow with
repetitive overlaps; this change does not bound the entire merge pipeline.

## Reproduction And Tests

The optional benchmark harness measures native processes and hashes their output
before removing its own generated edge files:

```sh
python3 tests/benchmark_merge_dedup.py --core build/megahit_core \
  --work-dir /tmp/merge-dedup-benchmark --long-prefix READ_COUNT_PREFIX \
  --contig SHORT_CONTIG_FASTA --long-contig LONG_CONTIG_FASTA \
  --short-k 21 --long-k 29 --threads 1,8 --repeats 3 --sort-memory 128
```

Per-run timings, peak RSS, temporary traffic, and output hashes are saved in
`results.json`. Python and `/usr/bin/time` are only benchmark-harness dependencies.

`cmake --build build --target merge_dedup_test` runs the exact sequence oracle,
all packed widths through k255, hash/FM overlap indexes, duplicate/palindromic
inputs, empty inputs, forced multi-pass spills, ordinal restoration, and
temporary-write failure cleanup. Native integration tests cover selection,
saved options, Python-free execution, and resume.
