# Contig-End Mapping and Consensus

The production workflow is now native C++ and does not require Python. See
[native executable and resume details](native-executable.md). The historical
Python driver remains in the source tree only for regression comparison;
its numbered checkpoints cannot be resumed by the native driver.

After descending reduction assembly ends at k21 by default, its consensus
local contigs, full contigs, bubbles, and additional contigs seed a fresh k21
SdBG. This seed graph does not reuse first-pass reduced or read-counted edges. The flanking
pass then assembles the full configured k-list, by default
21,31,41,61,81,101,141. First-pass files remain unchanged.

At the handoff and every subsequent flanking transition, the consensus stage runs:

1. Extract contig ends to `OUT_DIR/tmp/k21/k21.contig_ends.fa` for the default handoff.
2. Map merged and unmerged FLASH gzip FASTA exports listed in `assembly_inputs.json` to that file with Minimap2.
3. Build outward consensus extensions in native C++.
4. Publish extended records as `intermediate_contigs/k21.local.fa` and its info file.
5. At the handoff, build the flanking seed graph from the smallest reduction k's sequences. Thereafter,
   extract iterative edges and build the next SdBG with the new local-contig
   input, following the configured flanking list through k141. With
   `--flanking-consensus-only`, no read-derived iterative edges are added.

The handoff uses the paths shown above. The second pass uses
`intermediate_contigs/flanking/kN.*` and `tmp/flanking/kN/` to avoid overwriting
the first pass. Final contigs are merged from the second pass only. Exported
contig ends, consensus files, and alignments survive normal scratch cleanup
in both passes.

FLASH preprocessing is enabled by default and runs before binary-library
construction. Merged sequences and both unmerged mates are used as single-end
assembly inputs, and the same exports supply consensus. `--no-merge-pairs`
disables this preprocessing and FLASH-based consensus; without these inputs the
local FASTA is empty and its info file is `0 0`. There is no automatic fallback to
IDBA. `--no-local` disables local-contig publication and graph input entirely.
The old core `local` subcommand remains available for standalone use but is
not called by the launcher. Custom lists start flanking at the first flanking k,
which can differ from the smallest reduction k used as the source. Workflows
without a reduction list retain their increasing-k path without a handoff.

## Files

```text
tmp/k21/k21.contig_ends.fa
tmp/k21/k21.contig_ends.extended.fa
tmp/k21/k21.contig_ends.extended.fa.info
tmp/k21/flank_mapping/library.N.sam
tmp/k21/flank_mapping/mapping.json
tmp/k21/flank_mapping/consensus.tsv
tmp/k21/flank_mapping/consensus.json
intermediate_contigs/k21.local.fa
intermediate_contigs/k21.local.fa.info
```

The default handoff source is `intermediate_contigs/k21.contigs.fa`. Contigs >=580 bases produce
`original_name|start` and `original_name|end` fragments of 290 bases. Shorter
contigs retain the entire original sequence and header. Extracted fragments
remain in the original contig orientation. All output above remains under the
output directory even with a custom scratch directory and survives normal
temporary-file cleanup.

The extended FASTA contains only records with at least one added base. Original
flank sequence is retained inside each emitted record; records with no extension
are omitted. Headers begin with `flag=0 multi=N` for graph-reader compatibility,
where N is the number of reads contributing outward votes, capped at 65,535.
They also report uncapped `support_reads`, `reference_length`, `left_extension`,
`right_extension`, and the new `len`. This is one multiplicity per record, not
per-base pileup depth. The TSV still covers all reference fragments
for diagnostics, including those omitted from the FASTA.

`k21.local.fa` and its `.info` file are copies of the extended FASTA and its
`.info`. Each info file stores `<record_count> <sum_of_sequence_lengths>` followed
by a newline. Zero extensions produce an empty FASTA and `0 0` metadata. The
next graph builder applies its usual minimum sequence length of next_k + 1.

## Consensus Algorithm

The `endconsensus` core subcommand streams SAM and uses the project's two-bit
`SeqPackage` representation for read sequences, with native 64-bit A/C/G/T
counters for each overhang position. SAM text is parsed once; no Python
per-base pileup or external consensus package is used. References and FASTA
output remain strings to preserve their original sequence, including any Ns.

The native workflow forwards its CPU count to `endconsensus -t N`. Standalone
calls retain a one-worker default when `-t` is omitted. `--threads N` and
`--num-cpu-threads N` are also accepted before the positional arguments.

SAM files are split into byte ranges at query-group boundaries, so one large
SAM can use multiple workers without splitting a read's deduplication state.
Each worker uses its own input stream and packed-read scratch space. Shared
vote arrays/counters use 4,096 sharded locks rather than one complete pileup
per worker. Final consensus strings are computed in bounded blocks of 1,024
references, then FASTA, TSV, and info records are written in reference order.
Thread count does not change the consensus rules or output bytes. Worker parse
errors are collected before output files are opened; error locations in later
ranges use byte offsets. Parallel input files must support seeking.

On the saved E. coli k21 mappings, median consensus time over three repeats was
0.391 s before and 0.103 s with eight workers (3.78x). A single-SAM synthetic
workload with one million query groups took 2.169 s before and 0.375 s with
eight workers (5.79x). These are consensus-only, cached-file measurements, not
whole-assembly speedups. All FASTA, TSV, and `.info` hashes matched exactly.
Reproduce with `tests/benchmark_consensus.py`; performance depends on I/O,
reference distribution, and contention at highly covered ends.

Only mapped primary records vote. Secondary and supplementary alignments are
excluded to avoid repeatedly voting from alternative alignments of a read.
Reverse-strand SAM SEQ is already reverse-complemented by the aligner; it must
not be reverse-complemented a second time.

Leading and trailing soft clips are projected onto reference coordinates from
their aligned anchors. Only bases projected outside the reference contribute.
CIGAR M, =, and X consume both coordinates; I consumes query coordinates only;
D/N consumes reference coordinates only. H/P add no observed bases. Internal
insertions/deletions therefore do not shift the trailing-clip coordinate
incorrectly, but they are not themselves assembled into an indel consensus.

`|start` fragments extend to the left only; `|end` fragments extend to the right
only. Intact short contigs can extend on both sides. This avoids duplicating
known internal contig sequence at the artificial cut edges of the 290-base
fragments.

From each physical boundary outward, the base with the largest count wins.
Five A votes versus four G votes yields A. A unique plurality is sufficient;
no majority fraction or quality weighting is imposed. N and other ambiguous
read bases do not vote. Extension stops at the first uncovered position or
tie for the largest count. Bases beyond that stop are not appended.

## Scope and Limitations

This is a naive, single-pass overhang consensus, not a full variation-aware
assembly or polishing tool. One supporting read is sufficient in this initial
implementation. MAPQ is not filtered, base qualities are not weighted, and
misaligned or repetitive reads can produce incorrect extensions. There is no
iterative remapping beyond the reach of the currently aligned read overhangs.

Extended fragments are not spliced back into full original contigs. Instead,
they replace the standard local-contig contribution to the next SdBG. The
exact iterative-edge extractor still uses original reads and contig/bubble
sequences; `reads.lib.bin` remains unchanged.

Mapping and consensus have independent completion records. Input/reference
changes invalidate reuse; unchanged extraction preserves its timestamp so
mapping reuse remains possible. This does not fix the separately documented
pre-existing full-workflow reduction-stage resume issue. Use a fresh output
directory to apply the new local algorithm to runs completed with legacy local
assembly checkpoints. The extended-only algorithm version invalidates older
unfiltered consensus completion records when consensus is executed again.

## Standalone Use and Tests

```sh
megahit_core endconsensus ends.fa extended.fa consensus.tsv library.0.sam library.1.sam
cmake --build build --target end_consensus_test
```

The command pools votes across all supplied SAM files and writes
`extended.fa.info` automatically. The launcher publishes all output artifacts
after successful processing.
