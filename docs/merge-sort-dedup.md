# Exact Sort-Based Merge Deduplication

The default merge backend replaces both merge Bloom filters with a native,
bounded-record-buffer external sort:

```sh
megahit --merge-dedup sort --merge-sort-memory 128 [normal input/output options]
```

`--merge-sort-memory` is the total record-buffer budget in MiB, shared by sorting
workers, not a per-worker budget. Its default is 128 MiB and its minimum is
0.0625 MiB. The budget is capped at one quarter of the smaller visible available
memory and workflow `-m` budget. Indexes, overlap matches, I/O buffers, merge
cursors, and runtime overhead are additional; this is not an enforced process
RSS limit. The overlap index remains independently selectable with
`--merge-index hash|fm`.

Bloom remains selectable with `--merge-dedup bloom`.
The selection is saved in `options.json`; resume uses saved options and existing
checkpoints. Older native runs without a saved selection retain historical Bloom
deduplication. Use a fresh output directory to compare backends. The option applies
only to assembly-edge merging, not to the separate k-mer reduction Bloom filter.
The installed native executable still requires no Python, GNU sort, or other
external sorting program.

## Method

1. Copy every read-counted base edge unchanged to its original output shard.
2. Collect its canonical packed key with ordinal zero into per-worker buffers.
3. Generate extension candidates using the existing exact-overlap method. Store
   canonical keys with monotonically increasing generation ordinals.
4. Sort full buffers into binary runs using `std::sort`, discarding equal-key
   duplicates within each run. Flush remaining worker buffers in parallel.
5. Merge runs by key and ordinal. Base ordinal zero takes precedence, suppressing
   any candidate already present in the input. Otherwise retain the earliest
   candidate for each key. Multi-pass merging opens at most 32 input runs at once.
6. Sort accepted candidates by their retained ordinals and emit them in original
   generation order, with the existing multiplicity 50 and round-robin shard
   assignment.

Records store only active 64-bit sequence words plus one ordinal: 16 bytes at
k29 and 24 bytes at k39/k59. Neither merge Bloom filter is allocated. Base input
multiplicities, orientations, and even pre-existing duplicate base records are
preserved; this backend makes newly added edges unique, rather than rewriting
the base library or aggregating its counts.

The external sorter is in `bloom_filter/exact_merge_dedup.cpp`. Temporary run
directories are siblings of the merged edge prefix, for example
`tmp/k29/merged_edges/29.sort-XXXXXX`. They are removed on success and handled
exceptions. Abrupt termination can leave a workspace behind. Run read/write
failures are checked, and incomplete work does not publish new edge metadata.
Per-candidate provenance text reports are not supported with sort mode; the core
requires `--human_readable 0`, as the native workflow already supplies.

Sorting is exact, so Bloom false positives cannot discard candidates. Byte-for-
byte equality with Bloom is expected when its membership answers are correct,
but is not guaranteed on arbitrary data because Bloom can have false positives.

## Benchmarks

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
