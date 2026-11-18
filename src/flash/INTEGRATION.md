# FLASH Integration

Upstream: https://github.com/ebiggers/flash
Revision: 304e456baea2696e85b866752dd10984668f603d (FLASH v1.2.11)
License: GPL-3.0-or-later; see COPYING.

The upstream C sources are unmodified. CMake renames the FLASH entry point
to `main_flash` and links it into `megahit_core`. Run
`megahit_core flash --help` for the original FLASH command-line interface.
No separate FLASH executable or network access is required at runtime.

The native workflow runs FLASH before binary-library construction, enabled by
default for paired FASTQ (`--merge-pairs` is accepted explicitly). Successfully
merged sequences and both unmerged mate outputs enter `reads.lib.bin` as
single-end libraries. Original mates of merged pairs are not added again.
Use `--no-merge-pairs` for the previous original-mate workflow. The binary format
is unchanged. The historical Python driver remains a regression reference;
production workflow behavior is documented in [native details](../../docs/native-executable.md).
Separate and interleaved paired FASTQ
inputs are supported, including gzip and bzip2 compression. Paired FASTA
inputs are not supported in this mode because FLASH requires qualities.

Parameters exposed by the launcher are `--flash-min-overlap` (10),
`--flash-max-overlap` (65), `--flash-max-mismatch-density` (0.25), and
`--flash-phred-offset` (33). The maximum overlap affects scoring, not a
hard upper limit on accepted overlaps. FLASH uses `-t` from the launcher;
with multiple threads output order may vary, but mates remain associated.

Each paired library produces `OUT_DIR/tmp/flash/library.N.merged.fa.gz`,
containing successfully merged sequences, plus `library.N.unmerged.1.fa.gz` and
`library.N.unmerged.2.fa.gz`. All are gzip FASTA. Nonempty exports are listed in
`assembly_inputs.json` and used by both binary-library creation and consensus.
Conversion streams validated FASTA records to bundled `pigz -1 -c -n -p NUM_THREADS`
when the optional compressor was built and `-t` is greater than one. Compression is parallel;
FASTQ parsing/FASTA formatting remains a single streaming producer, with
buffered writes and no intermediate uncompressed FASTA. Without pigz, or
with `-t 1`, native C zlib gzip at level 1 is used. Pigz source is vendored in
`src/pigz` and built/embedded by default. If that optional build fails, MEGAHIT
still builds and uses zlib gzip; system pigz is not required or searched.
Use `-DBUNDLED_PIGZ=OFF` to omit it. See [build details](../pigz/INTEGRATION.md).
Both paths preserve identical decompressed
sequence content; gzip bytes and file sizes may differ. Compressor selection,
thread count, and conversion elapsed time are logged. A failed compressor is
reaped and its partial output removed; the export and completion manifest
are published only after successful conversion.
Libraries are numbered from zero, with separate-file pairs first, followed
by interleaved inputs. A library with no mergeable pairs produces a valid
empty gzip file. Single-end inputs are not included in the merged export.

The output location remains under `OUT_DIR/tmp`, even with `--tmp-dir`.
FASTA exports and the atomic `merged.json`/`assembly_inputs.json` manifests survive normal
temporary-file cleanup; intermediate FLASH FASTQ files are removed.
Native resume reuses completed exports or regenerates missing exports. A new
preprocessing checkpoint invalidates older binary libraries and downstream stages
on incomplete FLASH-enabled runs, so original mates cannot be silently reused.

With merged export enabled, the [Minimap2 stage](../minimap2/INTEGRATION.md)
maps both merged and unmerged exports against extracted 290-base contig ends at every flanking
transition. Mapping writes per-current-k SAM
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
