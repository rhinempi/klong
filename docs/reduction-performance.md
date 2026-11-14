# Edge Reduction Performance

Reduction still removes a short edge only when both its forward sequence and
reverse complement hit the longer-edge prefix Bloom filter. Longer edges insert
their forward prefix and reverse-complement suffix. Multiplicities and binary
edge encoding are unchanged. Normal execution allocates no provenance hash map
and writes no per-edge text reports.

## Implementation

The optimizations were implemented in this order:

1. Retained edges stream into bounded 1-MiB output buffers, eliminating the
   vector-of-vectors and individual retained-edge heap allocations.
2. Input records are read in 1-MiB blocks through an inlined templated callback.
   Missing shards, incomplete records, and failed output writes are errors.
3. Prefix extraction and reverse complementation operate on packed words.
   Reduction hashes only the active 64-bit words: one for k21/k29, two for k39.
   A failed forward lookup skips the reverse lookup when diagnostics are off.
4. Bloom sizing targets an estimated false-positive rate of `1e-15`. It sizes
   for up to 12 preferred probes (bounded by `--max_hashes`), then selects the
   smallest probe count meeting the target in the allocated space. Power-of-two
   allocations use bit masks instead of division. If rounding up would exceed
   the cap, the full aligned cap is available rather than discarding up to half
   the budget. Initialization is parallel. The flat Bloom layout is retained;
   no blocked-filter approximation was introduced.
5. Reduction receives the workflow's `-m` budget as `--host_mem`. Allocation also
   considers host available memory and visible cgroup v1/v2 limits. It budgets
   2 MiB per worker for input/output buffers and leaves a reserve of up to
   `--reserve_gb`, limited to one eighth of the available budget. This leaves
   room for metadata/runtime overhead but is not an OS-enforced RSS limit.
6. Large input files are divided into record-aligned chunks, targeting roughly
   four chunks per worker with a minimum chunk size of 1 MiB. A single large
   shard can therefore use multiple workers. Output streams remain private to
   their worker. Unsorted record order and shard assignment can change; edge
   multisets, not shard bytes, are the appropriate comparison.

The Bloom filter remains probabilistic. Its size, probe count, active word count,
estimated false-positive rate, and phase timings are logged. If the memory or
hash cap prevents the target rate, a warning is emitted and the best available
probe count is used. The target is a standard Bloom-model estimate, not an exact
membership guarantee. The new hash/probe settings need not reproduce historical
false positives.

Standalone reducer options include `--host_mem BYTES`, `--reduce_fpr RATE`,
`--max_bloom_gb GIB`, `--reserve_gb GIB`, and `--max_hashes N`. The native workflow
automatically passes its memory budget and uses the default FPR target. Neither
production builds nor runtime require Python.

## Validation

`tests/test_reduce.py` compares reduction against an independent sequence oracle
and optionally an older executable (`MEGAHIT_TEST_BASELINE_CORE`). It covers
both orientations, multiplicities, duplicates, empty inputs, packed boundaries
through k255, sorted/unsorted metadata, input/output buffer boundaries, multiple
workers on one large shard, capped filters, and malformed inputs. Run it with:

```sh
cmake --build build --target reduction_test
```

Python is only the optional test harness. Native integration tests check that the
workflow passes its memory budget into reduction and still skips preliminary
read-counted SdBG construction.

## E. coli Benchmark

Three alternating before/after runs per transition used the same retained
E. coli binary edges, eight workers, warm filesystem caches, and no diagnostic
reports. The optimized process received a 4,000,000,000-byte budget; the old
process used its original 60-GiB Bloom cap. Both fit within 4 GB on this dataset.
Numbers below are medians of complete process elapsed time and peak RSS.

| Reduction | Before | After | Speedup | Before RSS | After RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| k21 using k29 | 2.53 s | 1.18 s | 2.14x | 566.9 MiB | 274.1 MiB |
| k29 using k39 | 2.69 s | 1.18 s | 2.28x | 566.4 MiB | 270.8 MiB |
| k39 using k59 | 2.65 s | 1.17 s | 2.26x | 567.8 MiB | 270.9 MiB |

Every run retained the same edge multiset as the baseline: 6,526, 12,594, and
32,860 edges respectively. The optimized filters used 256 MiB and 12 probes,
with estimated FPRs around `2.4e-16`. These results measure reduction only;
speedups on larger metagenomes depend on storage, memory bandwidth, thread
count, edge distribution, and whether the Bloom memory cap is reached.
