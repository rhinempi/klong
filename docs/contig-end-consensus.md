# Contig-End Mapping and Consensus

At every flanking transition (starting at k59 for the default k-list), the
replacement consensus local stage runs:

1. Extract contig ends to `OUT_DIR/tmp/k59/k59.contig_ends.fa`.
2. Map merged gzip FASTA libraries to that file with Minimap2.
3. Build outward consensus extensions in native C++.
4. Publish extended records as `intermediate_contigs/k59.local.fa` and its info file.
5. Extract iterative edges, build the k79 SdBG with the new local-contig input,
   and assemble it. Repeat the same process from k79 to k99 and later phases.

Mapping and consensus require `--merge-pairs`; without merged reads, the local
FASTA is empty and its info file is `0 0`. There is no automatic fallback to
IDBA. `--no-local` disables local-contig publication and graph input entirely.
The old core `local` subcommand remains available for standalone use but is
not called by the launcher. The optional reduction-handoff export still runs
without a subsequent transition when `--merge-pairs` is supplied.
Custom k-lists and workflows starting above k60 use their actual current k.

## Files

```text
tmp/k59/k59.contig_ends.fa
tmp/k59/k59.contig_ends.extended.fa
tmp/k59/k59.contig_ends.extended.fa.info
tmp/k59/flank_mapping/library.N.sam
tmp/k59/flank_mapping/mapping.json
tmp/k59/flank_mapping/consensus.tsv
tmp/k59/flank_mapping/consensus.json
intermediate_contigs/k59.local.fa
intermediate_contigs/k59.local.fa.info
```

The source is `intermediate_contigs/k59.contigs.fa`. Contigs >=580 bases produce
`original_name|start` and `original_name|end` fragments of 290 bases. Shorter
contigs retain the entire original sequence and header. Extracted fragments
remain in the original contig orientation. All output above remains under the
output directory even with a custom scratch directory and survives normal
temporary-file cleanup.

The extended FASTA contains only records with at least one added base. Original
flank sequence is retained inside each emitted record; records with no extension
are omitted. Headers begin with `flag=0 multi=1` for graph-reader compatibility
and also report `reference_length`, `left_extension`, `right_extension`, and
the new `len`. Multiplicity 1 follows the previous local-contig convention;
it is not a measured pileup depth. The TSV still covers all reference fragments
for diagnostics, including those omitted from the FASTA.

`k59.local.fa` and its `.info` file are copies of the extended FASTA and its
`.info`. Each info file stores `<record_count> <sum_of_sequence_lengths>` followed
by a newline. Zero extensions produce an empty FASTA and `0 0` metadata. The
next graph builder applies its usual minimum sequence length of next_k + 1.

## Consensus Algorithm

The `endconsensus` core subcommand streams SAM and uses the project's two-bit
`SeqPackage` representation for read sequences, with native 64-bit A/C/G/T
counters for each overhang position. SAM text is parsed once; no Python
per-base pileup or external consensus package is used. References and FASTA
output remain strings to preserve their original sequence, including any Ns.

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
