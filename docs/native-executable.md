# Native Executable

The installed `megahit` is a compiled Linux ELF executable. Its workflow
driver is C++11 and its assembly core is linked directly into the same file.
It does not execute or link Python, extract a Python launcher, or require a
separate installed `megahit_core`. Build-time Python is not required either.
The optional Python test harness and historical `src/megahit` driver remain
in the source tree for regression comparisons, not production execution.

## Dependencies

Build dependencies are a C/C++11 compiler, CMake 3.5+, OpenMP, zlib development
headers, and bzip2 development headers. On Ubuntu/Debian:

```sh
sudo apt-get install -y git build-essential cmake zlib1g-dev libbz2-dev
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -- -j4
./megahit --version
./megahit --native-info
```

Normal builds are dynamically linked against standard C/C++ runtime, OpenMP,
zlib, and bzip2 libraries. "One executable" does not mean a completely static
binary independent of the OS libraries. `STATIC_BUILD=ON` retains the project's
static-link option where static dependency libraries are available.

Only `megahit` is installed. `megahit_core` is available as an explicit
development/test build target but the production workflow invokes its own executable with `--core`.
Native subcommands can also be run as `megahit minimap2 ...`, `megahit flash ...`,
or `megahit endconsensus ...`.

## Compression

Pigz's vendored C source is compiled in an isolated optional CMake sub-build.
The resulting binary is embedded as bytes in the native executable and unpacked
to an owned temporary directory only when parallel compression is requested.
That temporary directory is removed on normal completion and handled failures.
The assembly core is linked, not extracted. No system pigz or PATH lookup is used.

If pigz cannot build, `BUNDLED_PIGZ=OFF` is set, unpacking is unavailable, or only
one thread is requested, gzip FASTA is written using native zlib at level 1.
There is no Python gzip fallback. Failed compression removes the partial export.
Gzip/bzip2 input decompression also uses native libraries, streamed through owned
FIFOs when constructing the binary read library. FASTA formatting is native.

## Read Preprocessing

Paired FASTQ input runs through FLASH before `reads.lib.bin` is built. This is
enabled by default (`--merge-pairs` remains accepted). Each pair contributes either
one merged sequence, or its two unmerged mates as separate single-end sequences.
The original mates of a successfully merged pair are not added again. Original
`-r` single-end libraries are also retained in the binary library.

Each paired library exports gzip FASTA under `OUT_DIR/tmp/flash/`:

- `library.N.merged.fa.gz`
- `library.N.unmerged.1.fa.gz`
- `library.N.unmerged.2.fa.gz`

`assembly_inputs.json` lists nonempty exports used both for binary-library
construction and for consensus mapping at every flanking transition. Empty exports
are still valid gzip files, but are omitted from that list. `merged.json` separately
lists merged exports for provenance, including empty ones. Consensus does not add
the original `-r` libraries. All FLASH-derived binary libraries are marked single-end.
These exports/manifests survive normal cleanup, including with custom scratch.

Separate and interleaved paired FASTQ, gzip, and bzip2 input are supported.
`--no-merge-pairs` retains the original mate-based library behavior and disables
FLASH-based consensus. Use it for paired FASTA, which FLASH does not support.
The built-in `--test` uses bundled FASTA libraries and skips FLASH.

Incomplete older FLASH-enabled runs rebuild their exports, binary library, and
downstream stages under the versioned preprocessing checkpoint. Completed runs
remain completed on `--continue`; use a fresh output directory to rebuild them.

## Workflow And Resume

Reduction/merge defaults to `21,29,39,59`; consensus flanking defaults to
`21,31,41,61,81,101,141`. Set them independently:

```sh
megahit -1 R1.fq.gz -2 R2.fq.gz --merge-pairs -o output \
  --reduction-k-list 21,29,39,59 \
  --flanking-k-list 21,31,41,61,81,101,141
```

Each list must be strictly increasing, with odd values in [15,255]. Adjacent
reduction steps may be at most 28; flanking steps may be at most 40. Explicit
reduction lists are not restricted to k<=60. The last reduction assembly and its
consensus extensions seed a fresh graph at the first flanking k, which may differ
from the first reduction k. The same 290-base end extraction, Minimap2
parameters, extended-only consensus local files, and binary edge formats are used.
Final contigs are filtered and merged in C++ without shell pipelines.

For compatibility, `--k-list`, `--k-min/--k-max/--k-step`, and presets retain their
shared-list behavior: reduction uses entries <=60 and flanking uses the full list.
A phase-specific list overrides that phase regardless of CLI option order.
Explicit phase lists and the new defaults are not shortened based on read length.
Both lists are saved in `options.json`; continuation restores them. Older native
runs without these fields retain their historical shared-list behavior.

`--flanking-consensus-only` skips read-derived iterative `.edges.*` generation
throughout the flanking pass. Each new SdBG instead uses retained contigs,
bubble/additional contig sequences, and consensus `.local.fa` sequences. Stale
iterative edge files are ignored. Retained contigs can carry an iteration forward
even without new consensus extensions. Reduction and the initial seeding graph
are unchanged. The option is off by default and saved for native resume; use a
fresh output directory to change an existing run's mode. Consensus uses both
merged and unmerged FLASH exports. Without paired libraries, with
`--no-merge-pairs`, or with `--no-local`, this mode can only rebuild from retained
contig sequences.

Contigs shorter than 580 bp are retained whole and may extend on both sides.
Contigs of 580 bp or longer produce separate 290 bp `|start` and `|end` references,
each eligible only for outward extension.

Consensus mappings use `-x sr -a -Y --secondary=no -t N -I 55G`.
Only primary alignments contribute outward soft-clipped bases; unmapped,
secondary, and supplementary alignments are excluded. No additional MAPQ or
minimum read-support threshold is applied. Ambiguous reads can still contribute
at the primary placement selected by Minimap2, including MAPQ 0 placements.

Each extended FASTA record has `multi=N`, where N counts reads contributing valid
outward A/C/G/T votes for that reference end, including minority votes. Reads
merely aligning inside the reference are not counted. This is one multiplicity
for the whole extended record, not per-base extension coverage. It is capped at
65,535 for binary edge storage; `support_reads` preserves the uncapped count.
`consensus.tsv` retains primary/overhang counts; the compatibility column
`secondary_reads` is now always zero. `primary_reads` counts alignment records;
`overhang_reads` is the per-query-group deduplicated support count used for multiplicity.

The consensus engine uses the workflow CPU count for query-boundary SAM ranges
and bounded parallel consensus calling; ordered output is byte-identical across
thread counts. Standalone syntax is `megahit endconsensus -t N ends.fa extended.fa
stats.tsv input.sam ...` (omit `-t` for one worker). Cached consensus remains valid
because threading does not change the algorithm or output format. See
[consensus details](contig-end-consensus.md) for memory behavior and measurements.

Alignments must be grouped by query name, as in the workflow's unsorted Minimap2
SAM output. A query contributes at most once per reference end within its group,
using the first primary placement with an eligible outward vote. Normal mapping
therefore uses a read at only its selected primary reference end. This uses only
per-query bookkeeping. Each SAM file is a
separate read library; query names should be unique within a library. Coordinate-
sorted or arbitrarily interleaved SAM must first be grouped by query name.

The versioned consensus checkpoint forces incomplete older runs to regenerate
consensus and downstream stages; changed mapping parameters invalidate cached
secondary-enabled SAM. Already completed runs still remain completed on `--continue`;
use a fresh output directory to rebuild them with the new consensus behavior.

Look-ahead counting retains only binary `.edges.*` files; it does not build a
preliminary read-counted SdBG. Reduction and merging consume those edges directly.
SdBG construction occurs only for the reduced/merged edges that will be assembled
and for the subsequent flanking assembly graphs. Native checkpoints from earlier
versions remain usable; obsolete read-counted graph checkpoints are ignored.
Reduction uses bounded block I/O, packed-word operations, and chunk-based workers;
its Bloom allocation now receives the workflow memory budget. See
[Edge Reduction Performance](reduction-performance.md) for parameters and benchmarks.
Assembly skips unused depth inference, uses thread-private histograms, avoids
no-change refreshes, and checks changed-only output before sequence reconstruction.
See [Assembly Performance](assembly-performance.md) for details and comparison commands.
Intermediate reduction assemblies output only contigs longer than the next k
(for example, at least 30 bp for k21 -> k29). The merge extracts packed
(next-k+1)-mer windows directly from those contigs and combines them with the
read-counted next-k edges. It no longer writes short/long split FASTA files or
searches short-contig ends against the next-k edges; no overlap index is built.
The final reduction assembly retains its usual output lengths. Its input is the
last `merged_edges` set, which includes the preceding assembly's novel edges;
this set is not generally redundant with `read_count`.

`--merge-index hash|fm` remains accepted for compatibility, but has no effect on
these long-only merges. The standalone core merge still supports short-contig
matching; see [Contig-End FM Index](merge-fm-index.md) for that mode's benchmarks.
The default `--merge-dedup sort` replaces merge Bloom filtering with exact external sorting;
see [Exact Sort-Based Merge Deduplication](merge-sort-dedup.md) for memory/disk
tradeoffs, measurements, and reproduction commands.
Intermediate merges now retain duplicate additions without filtering or sorting.
After reduction, only surviving records are exactly deduplicated into
`reduced_edges`, retaining the maximum multiplicity for each canonical key.
This is not a sum of counts. The terminal reduction k has no later reduction,
so its merge still uses `--merge-dedup`: exact sorting by default suppresses
duplicate additions and preserves original base records and counts. Parallel
extraction, worker-local sorting and key-range merging are described in the
deduplication document above. Versioned reduction and merge checkpoints rebuild the
affected stages of incomplete older runs. Completed runs remain completed on
`--continue`; use a fresh output directory to rebuild them with this behavior.

Native runs write `options.json` with `workflow_driver: native-cpp-v1` and
`checkpoints.native.json`. Resume validates input paths, sizes, and modification
times, then reuses completed stages with existing outputs. A missing stage output
causes that stage and subsequent stages to rerun. Stage artifacts are retained
until successful final cleanup for resume safety; peak temporary disk usage can
therefore be higher than with the historical per-round cleanup. `--keep-tmp-files`
preserves all scratch; otherwise exported consensus/mapping/merged FASTA files
survive the final cleanup as before.

Legacy Python checkpoint numbers are not compatible. A native executable refuses
to resume those runs rather than silently interpreting their checkpoints. Use a
fresh output directory when changing drivers. Existing running processes are
unaffected by replacing the installed executable.

## Validation

`tests/test_native_workflow.py` copies only the ELF into an isolated directory
and runs it with empty PATH and invalid Python environment variables. It checks
paired/interleaved compressed input, merged/unmerged binary reads, local consensus
publication, standalone output, native gzip/bzip2 decompression including
concatenated streams, and native resume. Set `MEGAHIT_TEST_NATIVE_ZLIB` to also
exercise an executable built with pigz disabled. Python is used only to drive
these tests; the spawned workflows do not use Python.

The native JSON manifests use vendored nlohmann/json 3.11.3 (MIT license in
`src/vendor/nlohmann/LICENSE.MIT`), avoiding an additional installed JSON library.
