# AI Handoff: MEGAHIT Klong Workspace

Workspace: `/mnt/klong/megahit`

Purpose: concise operational state for a future AI/code session.

## Current Status

- Branch: `master`
- Remote push target used successfully: `git@github.com:rhinempi/klong.git`
- SSH push may require:
  ```bash
  eval "$(ssh-agent -s)"
  ssh-add ~/.ssh/reflexiv
  git push origin master
  ```
- Last pushed commit known in conversation:
  - `a891cbf Require both orientations for prefix bloom reduction`
- Current uncommitted source change at time of note:
  - `bloom_filter/main_prefix_bloom.cpp`
  - merge mode scans contigs as `fwd`, `rc`, `rev`, and `comp`
- Many untracked generated dirs exist. Do not stage them:
  - `build-embedded/`
  - `build-single/`
  - `megahit_core`
  - `src/__pycache__/`
  - `test-*`

## High-Level Pipeline Implemented

Original MEGAHIT iteration was modified. Intended loop:

```text
count k21
count k29
reduce k21 using k29 prefixes
assemble k21 from reduced k21
merge/extend k21 contigs into k29
count k39
reduce k29 using k39 prefixes
assemble k29 from reduced k29
merge/extend k29 contigs into k39
...
reduce k119 using k141
assemble k119 from reduced k119
merge/extend k119 contigs into k141
assemble final k141
```

Important output prefixes:

```text
tmp/k21/21.edges.*                         raw first-k edges
tmp/kXX/read_count/XX.edges.*              raw read-counted long edges
tmp/kXX/reduced_edges/XX.edges.*           reduced/non-removed short edges
tmp/kXX/assembly_graph/XX.sdbg.*           SDBG from reduced edges
tmp/kXX/merged_edges/XX.edges.*            long read edges plus extension-derived edges
intermediate_contigs/kXX.contigs.fa        intermediate contigs
```

## Key Files

Main pipeline:

- `src/megahit`

Embedded build:

- `CMakeLists.txt`
- `cmake/embed_megahit.py`

Prefix Bloom reducer and merge/extension tool:

- `bloom_filter/main_prefix_bloom.cpp`

Edge reader/writer and orientation-sensitive internals:

- `src/sequence/io/edge/edge_reader.h`
- `src/sequence/io/edge/edge_writer.h`
- `src/sorting/kmer_counter.cpp`
- `src/sequence/io/binary_reader.h`
- `src/sequence/sequence_package.h`
- `src/sorting/seq_to_sdbg.cpp`
- `src/assembly/unitig_graph.cpp`

## Critical Orientation Facts

Do not forget this.

MEGAHIT stores reads reversed in counting paths:

- `kmer_counter.cpp`: `is_reverse = true`
- `binary_reader.h`: `AppendReversedCompactSequence()`
- `sequence_package.h`: plain reverse, not reverse complement

Downstream graph/output code is designed around this:

- `seq_to_sdbg.cpp` processes stored orientation and reverse complement.
- `unitig_graph.cpp::VertexToDNAString()` reverses labels before FASTA output.

Consequences:

- Edge text sidecars are stored/internal orientation.
- Contig FASTA output is external/original orientation.
- Direct `grep` of edge text against FASTQ can fail even when edge is correct.
- For matching contigs to internal edges, consider `fwd`, `rc`, `rev`, and `comp`.

Validated small-test edge correctness:

```text
expected_unique        497
observed_unique        497
missing_from_observed  0
extra_in_observed      0
multiplicity_mismatches 0
```

This was computed by reconstructing canonical k21+1 edges from reversed reads in `/mnt/klong/test_run/small_test.fq.gz`.

## Reduce Mode Semantics

Function: `RunReduce()` in `bloom_filter/main_prefix_bloom.cpp`

Invocation shape:

```text
megahit_core prefixbf --mode reduce \
  --short_prefix ... \
  --long_prefix ... \
  --output_prefix ... \
  --short_k K \
  --long_k NEXT_K
```

Current intended behavior:

- Operates on binary `.edges.*` records.
- Does not use text edge files.
- Edge length is k+1.
- Bloom filter indexes longer-edge prefixes only:
  - `prefix(long_edge)`
  - `prefix(reverse_complement(long_edge))`
- Do not index suffixes in reduce mode.
- A short edge is removed only if both short-edge orientations hit:
  - `forward(short_edge)`
  - `reverse_complement(short_edge)`
- One-orientation hits are kept.

Important example that drove a fix:

```text
short edge: GTTTCTTTCGCGCTAGAAACTA
long edge:  GCCTAATCGTTTCTTTCGCGCTAGAAACTA
RC(long):   TAGTTTCTAGCGCGAAAGAAACGATTAGGC
```

The short edge must be removed because:

- forward short hits `prefix(GTTTCTTTCGCGCTAGAAACTAAAGTTGCA)`
- reverse-complement short hits `prefix(RC(GCCTAATCGTTTCTTTCGCGCTAGAAACTA))`

Reports:

```text
*.removed.txt:
  #edge multiplicity reason matched_long_edge

*.non_removed.txt:
  #edge multiplicity reason
  reasons: forward_only_hit, reverse_complement_only_hit, no_orientation_hit

*.reduction.summary.txt:
  removed_short_edges
  non_removed_short_edges
  non_removed_forward_only_hits
  non_removed_reverse_complement_only_hits
  non_removed_no_orientation_hits
```

Small test after prefix-only correction:

```text
removed_short_edges      467
non_removed_short_edges  30
bloom_insertions         964
representation           binary_short_edges_reduced_by_raw_long_edge_and_reverse_complement_prefixes_both_short_orientations_required
```

## Post-Assembly Contig Split

Function: `split_contigs_by_next_k()` in `src/megahit`

Normal intermediate assemblies now split `intermediate_contigs/kK.contigs.fa` by `next_k + 1`:

```text
kK.lt_kNEXT_plus1.contigs.fa  length < next_k + 1
kK.ge_kNEXT_plus1.contigs.fa  length >= next_k + 1
```

`merge_contig_edges_into_next_k()` uses only the short split file as `--contig`. If it is missing, it falls back to `kK.contigs.fa` with a warning. The split is intentionally skipped when `--debug-assembly` disables pruning.

## Merge/Extension Mode Semantics

Function: `RunMerge()` in `bloom_filter/main_prefix_bloom.cpp`

Invocation shape:

```text
megahit_core prefixbf --mode merge \
  --long_prefix tmp/kNEXT/read_count/NEXT \
  --output_prefix tmp/kNEXT/merged_edges/NEXT \
  --contig intermediate_contigs/kK.contigs.fa \
  --short_k K \
  --long_k NEXT
```

Current intended behavior:

- Copies base long read-counted edges to merged output.
- Builds exact endpoint indexes over long edges:
  - prefix index
  - suffix index
- Long edges are indexed in stored orientation and reverse complement.
- Overlap length is exactly `short_k`, not `short_k + 1`.
- It searches useful endpoint pairs only:
  - contig prefix to long suffix
  - contig suffix to long prefix
- On a match, extend and emit all new long-k+1 windows.
- Extension-derived multiplicity is `50`.
- Dedupe with canonical long-edge key.

Important local/uncommitted orientation fix:

- Contigs are output in external/original orientation.
- Long edge text is stored/internal orientation.
- Therefore merge now scans contigs as:
  - `fwd`
  - `rc`
  - `rev`
  - `comp`
- Summary field added:
  - `contig_orientation_modes fwd,rc,rev,comp`

Small test after merge orientation fix:

```text
k21 -> k29:
contigs_seen                 6
contig_orientations_scanned  24
overlap_matches              20
extended_windows_seen        108
extended_edges_added         18
output_edges_total           500
final output                 1 contigs, total 518 bp
```

## Debug Assembly Flag

Flag:

```text
--debug-assembly
```

Effect:

- Intermediate assemblies run with pruning and next-k length threshold disabled.
- Final k assembly remains normal.
- Debug short contigs are in `intermediate_contigs/kXX.contigs.fa`.
- `kXX.final.contigs.fa` may be empty and should not be used alone to judge debug output.

## Commands

Build:

```bash
cmake --build build-embedded -j 2
```

Whitespace check:

```bash
git diff --check
```

Small FASTQ run:

```bash
MEGAHIT_EMBED_CACHE=/tmp/megahit-embedded-test \
./build-embedded/megahit \
  -r /mnt/klong/test_run/small_test.fq.gz \
  -o /tmp/megahit-test \
  --force --keep-tmp-files -t 2 --debug-assembly
```

Inspect reduction:

```bash
sed -n '1,40p' /tmp/megahit-test/tmp/k21/reduced_edges/21.reduction.summary.txt
sed -n '1,30p' /tmp/megahit-test/tmp/k21/reduced_edges/21.removed.txt
sed -n '1,30p' /tmp/megahit-test/tmp/k21/reduced_edges/21.non_removed.txt
```

Inspect merge:

```bash
sed -n '1,40p' /tmp/megahit-test/tmp/k29/merged_edges/29.merge.summary.txt
sed -n '1,30p' /tmp/megahit-test/tmp/k29/merged_edges/29.contig_edges.txt
```

## Editing/Commit Notes

- Prefer editing only tracked source files.
- Generated test/build outputs are untracked and should not be committed.
- Current docs were created to help resume work:
  - `docs/klong_handoff_human.md`
  - `docs/klong_handoff_ai.md`
- `apply_patch` has failed in this environment with:
  ```text
  bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted
  ```
  Prior edits used small Python scripts with escalated shell commands as a workaround.

## Recent Commit History

```text
a891cbf Require both orientations for prefix bloom reduction
3be32ce Add exact overlap edge extension index
61bfb0b Add debug assembly flag
22732e8 Tune binary prefix bloom filter for runtime
1d333b2 Add binary prefix bloom filter tool
324fc9a Add read-counted debug outputs for iterative ks
8233ddd Write readable edge sidecar files
829a218 Embed megahit_core into megahit launcher
```

## Next Recommended Action

If resuming from this exact state:

1. Review `git status --short`.
2. Confirm whether the merge-orientation fix in `bloom_filter/main_prefix_bloom.cpp` should be committed.
3. Run:
   ```bash
   cmake --build build-embedded -j 2
   git diff --check
   ```
4. Run the small FASTQ debug test.
5. Commit only:
   - `bloom_filter/main_prefix_bloom.cpp`
   - `docs/klong_handoff_human.md`
   - `docs/klong_handoff_ai.md`
