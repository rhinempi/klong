# FLASH Integration

Upstream: https://github.com/ebiggers/flash
Revision: 304e456baea2696e85b866752dd10984668f603d (FLASH v1.2.11)
License: GPL-3.0-or-later; see COPYING.

The upstream C sources are unmodified. CMake renames the FLASH entry point
to `main_flash` and links it into `megahit_core`. Run
`megahit_core flash --help` for the original FLASH command-line interface.
No separate FLASH executable or network access is required at runtime.

The Python launcher's `merge_paired_reads()` function runs FLASH as a separate
export after binary library construction when `--merge-pairs` is supplied.
It never changes input library paths or adds merged reads to `reads.lib.bin`.
Assembly uses all original reads, with the original pairing and binary format.
Separate and interleaved paired FASTQ
inputs are supported, including gzip and bzip2 compression. Paired FASTA
inputs are not supported in this mode because FLASH requires qualities.

Parameters exposed by the launcher are `--flash-min-overlap` (10),
`--flash-max-overlap` (65), `--flash-max-mismatch-density` (0.25), and
`--flash-phred-offset` (33). The maximum overlap affects scoring, not a
hard upper limit on accepted overlaps. FLASH uses `-t` from the launcher;
with multiple threads output order may vary, but mates remain associated.

Each paired library produces `OUT_DIR/tmp/flash/library.N.merged.fa.gz`,
containing only successfully merged sequences in gzip-compressed FASTA.
Libraries are numbered from zero, with separate-file pairs first, followed
by interleaved inputs. A library with no mergeable pairs produces a valid
empty gzip file. Single-end inputs are not included in the merged export.

The output location remains under `OUT_DIR/tmp`, even with `--tmp-dir`.
Merged FASTA and the atomic `merged.json` completion manifest survive normal
temporary-file cleanup; intermediate FLASH FASTQ files are removed.
Resume reuses completed exports or regenerates missing exports independently
of assembly checkpoints. No assembly checkpoint numbers are changed.

With merged export enabled, the [Minimap2 stage](../minimap2/INTEGRATION.md)
maps these files against extracted 290-base contig ends at every flanking
transition. Mapping preserves original reads and writes per-current-k SAM
files, e.g. `OUT_DIR/tmp/k59/flank_mapping`. Outward consensus extensions replace
the former IDBA local-contig graph inputs. The same merged library is reused
across phases rather than merged again at each k.

Build and test:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
cmake --build build --target flash_test
build/megahit --merge-pairs -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o merged_out
```
