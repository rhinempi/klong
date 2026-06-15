# Klong MEGAHIT Modification Handoff

This document summarizes the MEGAHIT changes we made in this working tree, why they were made, and the details that are easy to mix up. It is meant to be readable by a human picking this project back up later.

## Current Goal

The project is no longer trying to run MEGAHIT's original iterative strategy unchanged. We changed the iteration logic so every k value can be read-counted, compared against the next longer k, reduced, assembled, and then used to augment the next k.

The intended high-level loop is:

1. Count read edges for k21.
2. Count read edges for k29.
3. Reduce k21 edges using k29 prefixes.
4. Assemble k21 from the reduced k21 edges.
5. Use the k21 assembly to extend/augment k29 edges.
6. Repeat: count k39, reduce k29 using k39, assemble reduced k29, extend/augment k39, and so on.
7. Final k141 is assembled normally from the merged/augmented k141 edge set.

The reduced edge sets are supposed to drive subsequent graph construction. The raw read-counted edge sets are preserved for debugging/checkpointing.

## Important Orientation Context

This is the easiest source of mistakes.

MEGAHIT internally reads FASTQ sequences in reversed order in some core paths. This is not reverse complement; it is plain string reversal. The original code is designed around that convention:

- `src/sorting/kmer_counter.cpp` calls `SequenceLibCollection::Read(..., is_reverse = true)`.
- `src/sequence/io/binary_reader.h` calls `AppendReversedCompactSequence()` when `reverse=true`.
- `src/sequence/sequence_package.h` implements that as a reverse of base order, not reverse complement.
- `src/sorting/seq_to_sdbg.cpp` processes both stored orientation and reverse-complement orientation when building the graph.
- `src/assembly/unitig_graph.cpp::VertexToDNAString()` reverses labels before writing contig FASTA.

This means:

- Binary/read-counted edges are in MEGAHIT's stored internal orientation.
- Human-readable edge text files also show stored internal orientation.
- Contig FASTA output is in external/original orientation.
- A text edge may not grep directly against the original FASTQ. It may match after reverse, complement, or reverse-complement depending on which internal orientation was selected.

We validated this with `/mnt/klong/test_run/small_test.fq.gz`: reconstructing canonical k21+1 edges from reversed reads matched the raw edge text exactly.

## Major Code Changes

### Embedded executable

Earlier work changed the build so the Python `src/megahit` launcher can be embedded/wrapped with `megahit_core`, producing a single executable from `build-embedded/megahit`.

Relevant files:

- `CMakeLists.txt`
- `cmake/embed_megahit.py`
- `src/megahit`

### Read-counted debug outputs

The pipeline now creates read-counted edges for every k, not just the first k. Human-readable sidecar files were added for edge/debug inspection.

Typical output layout:

- `tmp/k21/21.edges.*`: raw first-k edges
- `tmp/kXX/read_count/XX.edges.*`: raw read-counted edges for later k
- `tmp/kXX/reduced_edges/XX.edges.*`: non-removed/reduced short-k edges
- `tmp/kXX/assembly_graph/XX.sdbg.*`: graph built from reduced edges
- `tmp/kXX/merged_edges/XX.edges.*`: next-k read edges plus extension-derived edges

### Binary prefix Bloom reducer

Implemented in:

- `bloom_filter/main_prefix_bloom.cpp`

Invoked through:

- `megahit_core prefixbf --mode reduce`

Current reduction logic:

- Longer edges are binary edge records, not text files.
- The Bloom filter indexes only prefixes of the longer edge orientation set:
  - `prefix(long_edge)`
  - `prefix(reverse_complement(long_edge))`
- It does not index suffixes for the reduction step.
- A short edge is removed only if both orientations of the short edge hit:
  - `forward(short_edge)` hits
  - `reverse_complement(short_edge)` hits
- If only one orientation hits, the short edge is kept.

Reports:

- `*.removed.txt` now includes `matched_long_edge`.
- `*.non_removed.txt` reports why a kept edge was kept:
  - `forward_only_hit`
  - `reverse_complement_only_hit`
  - `no_orientation_hit`
- `*.reduction.summary.txt` includes counters for those categories.

Important correction made later: suffixes were mistakenly indexed for Bloom reduction for a short time. That was wrong. Suffix/prefix matching belongs to the later extension step, not reduction.

### Post-assembly contig split

After every normal-mode intermediate assembly, the pipeline splits `kXX.contigs.fa` using the next k edge length threshold `next_k + 1`:

- `kXX.lt_kNEXT_plus1.contigs.fa`: contigs shorter than `next_k + 1`
- `kXX.ge_kNEXT_plus1.contigs.fa`: contigs longer than or equal to `next_k + 1`

The merge/extension step now uses only the short split file. This split is skipped in `--debug-assembly` mode because debug mode intentionally keeps pruning disabled for inspection.

### Exact overlap extension/merge

Implemented in:

- `bloom_filter/main_prefix_bloom.cpp`

Invoked through:

- `megahit_core prefixbf --mode merge`

The extension logic:

- Starts from the next longer k read-counted edge set.
- Copies those base edges to the merged output.
- Builds exact-match indexes over longer edges.
- Uses overlap length equal to the shorter k, not k+1.
- Searches only useful endpoint directions:
  - contig prefix to long-edge suffix
  - contig suffix to long-edge prefix
- On a match, it extends the longer edge using the short contig overhang.
- It emits all new longer-k+1 windows from the extended sequence.
- Extension-derived multiplicity is currently `50`.

The index implementation is a sorted flat vector of `OverlapEntry`, not a suffix tree/trie. This was chosen because exact fixed-length endpoint matching is faster and simpler with packed keys plus sorting/binary search than with a general string tree.

### Merge orientation fix

This change is currently local in the working tree unless it has been committed after this note was written.

Problem:

- Longer edges are stored in MEGAHIT internal reversed orientation.
- Short contigs from assembly are output in original/external FASTA orientation.
- The merge code previously scanned contigs only as `fwd` and `rc`.
- That missed cases where the long-edge index needed `reverse(contig)` or `complement(contig)`.

Current local fix:

- Long edges are indexed as stored orientation plus reverse complement.
- Contigs are scanned in four orientations:
  - `fwd`
  - `rc`
  - `rev`
  - `comp`

Small test evidence:

- Before this fix, extension matches were absent or ineffective.
- After the fix on `/mnt/klong/test_run/small_test.fq.gz`, k21 to k29 merge reported:
  - `overlap_matches 20`
  - `extended_windows_seen 108`
  - `extended_edges_added 18`
  - final output: `1 contigs, total 518 bp`

## Debug Flag

The launcher has a debug assembly mode:

```bash
--debug-assembly
```

This disables pruning and the next-k output threshold for intermediate assemblies so short unitigs/contigs can be inspected. Final k assembly stays normal.

Debug output is written to normal intermediate files such as:

- `intermediate_contigs/k21.contigs.fa`

Do not look only at `k21.final.contigs.fa` when debugging short contigs; it may be empty by design.

## Common Mixups

1. `k21` edge files contain 22-mers.
   Edge files store k+1 edges, not k-mers.

2. Human-readable edge strings are internal orientation.
   They may not match FASTQ directly.

3. Bloom reduction uses longer prefixes only.
   Suffixes are not part of reduction.

4. Extension uses both prefix and suffix endpoint indexes.
   That is separate from Bloom reduction.

5. Reverse and reverse complement are not interchangeable here.
   MEGAHIT storage uses plain reversal; biological graph logic also uses reverse complements.

6. The current working tree may include many generated directories.
   Ignore untracked build/test outputs unless explicitly cleaning them.

## Useful Commands

Build:

```bash
cmake --build build-embedded -j 2
```

Run small debug test:

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
sed -n '1,20p' /tmp/megahit-test/tmp/k21/reduced_edges/21.removed.txt
sed -n '1,20p' /tmp/megahit-test/tmp/k21/reduced_edges/21.non_removed.txt
```

Inspect merge:

```bash
sed -n '1,40p' /tmp/megahit-test/tmp/k29/merged_edges/29.merge.summary.txt
sed -n '1,30p' /tmp/megahit-test/tmp/k29/merged_edges/29.contig_edges.txt
```

Push with SSH key if agent is not running:

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/reflexiv
git push origin master
```

## Recent Commit Context

Recent commits at the time this was written:

- `a891cbf Require both orientations for prefix bloom reduction`
- `3be32ce Add exact overlap edge extension index`
- `61bfb0b Add debug assembly flag`
- `22732e8 Tune binary prefix bloom filter for runtime`
- `1d333b2 Add binary prefix bloom filter tool`
- `324fc9a Add read-counted debug outputs for iterative ks`
- `8233ddd Write readable edge sidecar files`
- `829a218 Embed megahit_core into megahit launcher`

## Recommended Next Checks

Before committing further logic changes:

1. Run `cmake --build build-embedded -j 2`.
2. Run `git diff --check`.
3. Run the small FASTQ debug test.
4. Check `k21 -> k29` merge summary for `overlap_matches` and `extended_edges_added`.
5. Check that generated output directories are not staged.
