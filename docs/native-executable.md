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
FIFOs when constructing the original read library. FASTA formatting is native.

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
fresh output directory to change an existing run's mode. Use `--merge-pairs` to
provide FLASH-merged reads for consensus; without merged reads or with
`--no-local`, this mode can only rebuild from retained contig sequences.

Look-ahead counting retains only binary `.edges.*` files; it does not build a
preliminary read-counted SdBG. Reduction and merging consume those edges directly.
SdBG construction occurs only for the reduced/merged edges that will be assembled
and for the subsequent flanking assembly graphs. Native checkpoints from earlier
versions remain usable; obsolete read-counted graph checkpoints are ignored.
Reduction uses bounded block I/O, packed-word operations, and chunk-based workers;
its Bloom allocation now receives the workflow memory budget. See
[Edge Reduction Performance](reduction-performance.md) for parameters and benchmarks.
Contig-end matching also supports `--merge-index fm` for an exact BWT/FM-index
backend. Hash remains the measured default; see
[Contig-End FM Index](merge-fm-index.md) for parallelism and comparison results.
The default `--merge-dedup sort` replaces merge Bloom filtering with exact external sorting;
see [Exact Sort-Based Merge Deduplication](merge-sort-dedup.md) for memory/disk
tradeoffs, measurements, and reproduction commands.

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
paired/interleaved compressed input, original binary reads, local consensus
publication, standalone output, native gzip/bzip2 decompression including
concatenated streams, and native resume. Set `MEGAHIT_TEST_NATIVE_ZLIB` to also
exercise an executable built with pigz disabled. Python is used only to drive
these tests; the spawned workflows do not use Python.

The native JSON manifests use vendored nlohmann/json 3.11.3 (MIT license in
`src/vendor/nlohmann/LICENSE.MIT`), avoiding an additional installed JSON library.
