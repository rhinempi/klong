# Minimap2 Integration

Upstream: https://github.com/lh3/minimap2
Release: v2.31, version 2.31-r1302
Revision: 3c28777e7e2dcc90f825de1b9f17a89cca7d4452
License: MIT; see LICENSE.txt.

Sources are vendored alongside FLASH. CMake renames upstream's C entry point
to `main_minimap2` and links the aligner into `megahit_core`. The embedded
MEGAHIT executable includes it without requiring a separate executable on PATH.
The x86 build preserves upstream SSE2/SSE4.1 runtime dispatch; ARM uses the
bundled NEON compatibility header (line endings normalized to LF).

Run `megahit_core minimap2 --version` or `megahit_core minimap2 --help` for
the bundled aligner's standalone command interface.

When `--merge-pairs` is enabled, `map_merged_reads(cur_k)` runs inside the new
consensus local stage before every flanking transition: k59 for 59->79, k79 for
79->99, and so on. Workflows starting above k60 use the same local stage for
their transitions. The legacy IDBA local command is no longer invoked by the
launcher. `--no-local` disables publication and use of local-contig graph
inputs; the optional reduction-handoff export remains available separately.

The reference is `OUT_DIR/tmp/k59/k59.contig_ends.fa`, extracted from the final
k59 reduction-phase assembly before mapping. It contains 290-base end fragments
or intact original contigs shorter than 580 bases, not final flanking assembly
contigs or bubble sequences.
Each merged gzip FASTA library is mapped independently as single-end reads:

```sh
megahit_core minimap2 -x sr -a -Y -t NUM_CPU_THREADS -I 55G \
  -o OUT_DIR/tmp/k59/flank_mapping/library.N.sam.tmp \
  OUT_DIR/tmp/k59/k59.contig_ends.fa \
  OUT_DIR/tmp/flash/library.N.merged.fa.gz
```

Successful SAM files are atomically renamed to `library.N.sam`. All libraries
must finish before `mapping.json` is published. That file records version,
parameters, input signatures, and output paths for resume. Changed inputs or
missing outputs trigger regeneration without changing assembly checkpoints.
No-reference/empty-reference assemblies are reported; an empty reference is
skipped. Empty merged libraries produce header-only SAM files.

Alignment output always lives under `OUT_DIR/tmp/k59/flank_mapping`, including
when `--tmp-dir` points elsewhere, and survives normal temporary-file cleanup.
The alignments now feed the [native outward consensus step](../../docs/contig-end-consensus.md).
Neither mapping nor consensus replaces the original binary read library.
Consensus now supplies the standard `k*.local.fa` graph input, replacing the
old local-assembly sequences; exact iterative-edge extraction is unchanged.
The `55G` argument is Minimap2's
index partition size, not a hard process memory limit.

```sh
cmake --build build --target minimap2_test
build/megahit --merge-pairs -1 reads_1.fq.gz -2 reads_2.fq.gz -o out
```
