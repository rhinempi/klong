# Assembly Hot-Path Optimizations

These changes preserve graph-cleaning thresholds and contig-extraction rules.
They do not introduce a different traversal or simplify the biological model.

## Changes

- Skip automatic minimum-depth inference when depth pruning is disabled
  (`--prune_level 0`). Weak-link disconnection does not use that threshold.
- Build minimum-depth multiplicity histograms in thread-private dense arrays,
  then reduce counts into the original histogram. Multiplicities are 16-bit;
  private counters use approximately 512 KiB per configured worker on 64-bit
  systems. The existing threshold-inference algorithm is unchanged.
- Build contig-length statistics in thread-private ordered maps and merge counts
  after workers finish. Long and sparse length distributions do not require a
  dense array indexed by contig length. Summary definitions are unchanged.
- Skip graph refresh after bubble removal, weak-link disconnection, or global
  low-depth removal when no mutations were recorded. Tip-removal guards use
  removals from the current threshold, not the cumulative count across thresholds.
  Existing local-depth pruning already guards refresh.
- In changed-only contig output, test the changed flag before reconstructing
  the nucleotide sequence.

These routines apply their own mutations before returning. No-change guards rely
on that contract; callers must not leave pending deletion/disconnection flags for
a later unrelated cleaning routine to apply.

## Regression Tests

Preserve the pre-change development core, then build the new one:

```sh
cp build/megahit_core /tmp/megahit-assembly-before
cmake --build build --target megahit_core megahit -- -j4
MEGAHIT_TEST_CORE="$PWD/build/megahit_core" \
MEGAHIT_TEST_ASSEMBLY_BASELINE=/tmp/megahit-assembly-before \
  python3 -m unittest discover -s tests -p test_assembly_optimizations.py -v
```

The test constructs a graph containing branching variants, tips, isolated paths,
and a cycle. It compares pruning levels 0-3, inferred and explicit depth cutoffs,
final and intermediate output, and 1/4 worker configurations. Inferred cutoffs,
statistics, sequences, flags, coverage metadata, and `.info` files are checked.
Generated IDs/order are ignored; reverse complements and circular rotations are
normalized for comparison.

## Runtime Comparison

```sh
python3 tests/benchmark_assembly.py \
  --before /tmp/megahit-assembly-before --after build/megahit_core \
  --graph SAVED_GRAPH_PREFIX --work-dir comparison-output \
  --profile reduction --threads 8 --repeats 1
```

`--profile reduction` reproduces the saved Sample0 k21 reduction assembly options:
one cleaning round, no tip/bubble/depth pruning, and weak-link disconnection.
`--profile cleaning` enables five rounds, bubbles, inferred depth, depth pruning
level 2, and changed-only additional-contig output. Use an appropriate saved graph
for each profile; this is an assembly-only comparison, not end-to-end runtime.

The script records commands, binary hashes, wall times, GNU time maximum RSS,
summary statistics, and per-file output counts/hashes in `results.json`.
Output normalization and external sorting are outside the timed assembly command.
Sorting uses GNU sort with a 128 MiB memory setting. This is a benchmark dependency,
not a production build/runtime dependency. Repeat order alternates; caches are not
flushed. A single pair is a smoke measurement, not a controlled cold-cache result.

## Existing Circular Coverage Variation

`UnitigGraph`'s unchanged loop-construction code initializes depth with the start
edge's multiplicity, then counts that edge again while walking the cycle. With
unequal multiplicities, selecting a different starting edge changes the reported
average. Parallel loop discovery can select different starts across runs. The
weighted-cycle regression test reproduces multiple coverage annotations in the
pre-change binary while verifying that the circular sequence is unchanged.

The benchmark defaults to strict normalized comparison, including coverage.
`--allow-loop-coverage-variation` explicitly permits only circular coverage
annotations to differ. It still compares all sequences, flags, record/base counts,
non-circular coverage annotations, and `.info` files. Strict hashes and circular
coverage annotations remain in `results.json`; `matches_reference` remains false
when coverage differs. A separate `sequence_and_nonloop_metadata_match` records
the narrower comparison. This is not full metadata equivalence.

Use `--recheck` with the original arguments to audit already-generated outputs
without rerunning assembly or replacing the original runtime measurements. The
circular coverage behavior is documented here, not changed by these optimizations.

## Measured Results

Eight-thread, assembly-only comparisons on saved graphs:

| Graph/profile | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| Sample0 k21, 339.8 million SdBG entries, reduction profile; one pair | 369.31 s | 240.26 s | 1.54x |
| E. coli k59, cleaning profile; median of three pairs | 8.46 s | 3.17 s | 2.67x |

Sample0 peak RSS was effectively unchanged at 2,587 MiB. Its sequences, flags,
counts, and non-circular metadata matched; two circular coverage annotations
differed due to the existing issue above. E. coli matched strictly across all
six runs, including inferred depth 3.162, coverage annotations, and `.info` files.
Sample0 graph-load times differed (19.94 vs 7.86 s); excluding graph loading,
the measured time reduction was 33.48%. Caches were not flushed.
