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
reduction lists are not restricted to k<=60. The lists are supplied in ascending
order; the reduction phase assembles them in descending order. The resulting
smallest-k assembly and its consensus extensions seed a fresh graph at the first
flanking k, which may differ from the first reduction k. The same 290-base end extraction, Minimap2
parameters, extended-only consensus local files, and binary edge formats are used.
Final contigs are filtered and merged in C++ without shell pipelines.

For compatibility, `--k-list`, `--k-min/--k-max/--k-step`, and presets retain their
shared-list behavior: reduction uses entries <=60 and flanking uses the full list.
A phase-specific list overrides that phase regardless of CLI option order.
Explicit phase lists and the new defaults are not shortened based on read length.
Both lists are saved in `options.json`; continuation restores them. Older native
runs without these fields retain their historical shared-list behavior.

### Experimental End Anchors

End anchoring is off by default. Standard flanking builds each SdBG from
the complete prior contigs, not stitched anchors. This also applies with
`--flanking-consensus-only`: consensus still uses read mapping at real contig
ends, but those mapping windows do not replace the full graph inputs. To
switch an earlier anchored run to the standard workflow, start in a fresh
output directory and omit the end-anchor and `--anchor-*` options.

`--flanking-end-anchors` makes flanking extend or join saved contigs rather than
reassemble their full lengths. For example, to preserve the k31 reduction
assembly while trying extension at k21:

```sh
megahit -1 R1.fq.gz -2 R2.fq.gz -o anchored \
  --reduction-k-list 31,41,51 --flanking-end-anchors \
  --flanking-consensus-only --keep-tmp-files
```

Each flanking end anchor is 290 nt. Contigs longer than 580 nt supply one
stitched record with the first and last 290 nt. Contigs at or below 580 nt
supply their full sequence to the same SdBG and can extend through consensus
edges. Circular contigs are retained unchanged. Reads are still
mapped to real contig ends. The extension SdBG receives stitched anchors,
trimmed outward consensus, and, for later iterations when enabled, iterative
edges indexed from the stitched anchors. It does not receive previous bubble
or additional-contig sequences.

The synthetic FASTA and raw SdBG unitigs are scratch files under
`tmp/flanking/kN/`; `--keep-tmp-files` retains them. Flanking uses the raw
unitigs by default. The matcher selects one best full hit per anchor, on
either strand, favoring fewer substitutions and then longer unitigs.
Ties are resolved deterministically rather than rejected. Non-overlapping
hits on the same unitig are restored together: each anchor is replaced with
its full original in the matching orientation, and the intervening sequences
are copied once. A contig without a selected placement is retained unchanged.
An unselected stitched-anchor junction or an ambiguous base can still stop a
bridge. There is no independent validation of joins. The usual final length
filter still applies.
Restoration packs DNA at two bits per base. Exact 16-base seeds identify
candidates for full packed comparisons; OpenMP scans batches of candidates
with up to `-t` workers. The best hit is chosen across all batches and both
strands before the restored FASTA is written in source order.
Repetitive anchors may be hard to recover intact at low k, so this mode can
preserve contiguity at the cost of accepting fewer extensions.

Add `--anchor-max-mismatches 0` to either anchor mode to report approximate
matches without changing restoration. Overlapping 12-base seeds find
candidates on both strands; packed comparisons count substitutions over the
ungapped overlap. `kN.anchor-matches.tsv` is written next to the restored
contigs. It reports the best full match and best fragment when found, including
the aligned length, substitutions, and unaligned bases separately. A fragment
requires at least 24 aligned bases, and reported matches have at least 90%
identity. A no-hit result can mean there was no shared seed, not that every
base mismatched. This search is also parallel by candidate sequence, but a
repetitive seed can generate many candidate offsets before deduplication.
It does not model indels or stitch across candidates.

Values `1` through `8` also allow restoration from a best full
placement within that substitution limit, provided both outer k-mers match
exactly. Partial hits never trigger restoration, and the original contig
sequence is always inserted unchanged.

For a gapped alternative, use `--flanking-end-anchors --anchor-min-identity 30`
with graph threading off. Exact matches still win at 100%. Otherwise, the
bundled Edlib aligner aligns the whole anchor, on both strands, to a substring
of each eligible raw unitig. The best edit-distance alignment from each search
is ranked by identity: identical bases divided by all alignment columns,
including substitutions and gaps. Restoration replaces the aligned **target**
span with the original contig, so insertions and deletions do not shift the
splice coordinates. `kN.anchor-alignments.tsv` reports identity, target span,
gaps, and whether the hit was restored. The mode is parallel by anchor but is
more expensive than exact matching. A 30% cutoff can accept unrelated
alignments and should be used experimentally. Do not combine this option with
`--anchor-max-mismatches` or `--anchor-graph-thread`.

`--anchor-graph-report` checks whether each complete anchor spells an exact
path in the cleaned SdBG, including across unitig boundaries. The assembler
looks up the first k-mer, then follows the surviving edge that matches each
next base. Queries run in parallel without building a separate unitig index.
The `kN.anchor-graph-paths.tsv` report beside the restored contigs shows the
matched prefix length and branch crossings. A full path does not establish
that an outward extension is unique. The report itself does not change restoration.
The report uses the graph snapshot that produces the primary `.contigs.fa`,
not the later `.addi.fa` for non-final iterations.

`--anchor-graph-thread` is off by default. It writes
both the graph path report and `kN.anchor-threaded.fa`. In flanking, only
stitched anchors are threaded. Restoration matches those against the threaded
FASTA, and matches whole contigs of 580 nt or less against the raw SdBG
unitigs. These shorter contigs are marked `raw_unitig` in the graph path
report. Short and long placements are evaluated separately; a complete
stitched-anchor match in a raw unitig also blocks copying through that anchor
during short-contig restoration.
In the threaded FASTA, several records can cover the same exact anchor path.
When one exact hit contains every other hit at the aligned anchor position,
restoration uses the anchor's own exported record and
does not count the other copies as competing placements. Conflicting paths
and repeated hits within one record remain ambiguous. Other records retain
their own anchors and are not removed from the FASTA.
Use `--no-anchor-graph-thread` to explicitly select raw-unitig restoration,
or `--anchor-graph-report` to keep only the TSV. Neither flag enables end-anchor
assembly by itself. A FASTA record is emitted only when the graph spells the
full anchor. Its sequence is spelled from the graph, with any left or right
flank extended through unique links until a branch, dead end, or edge reuse.
The FASTA retains synthetic stitches: only the restored contigs with their
original interiors reinserted are carried forward. In the descending
reduction phase, threaded paths are diagnostic only and restoration continues
to use the raw unitigs.

`--reduction-end-anchors` independently applies the same restoration to the
descending reduction pass. Ascending counting/reduction and the largest-k
assembly are unchanged. Each smaller k's SdBG uses the preceding larger-k
contig anchors and the target k's `reduced_edges`, without consensus or
iterative read edges. Anchor length is at least 100 nt and covers the maximum
of both configured k lists. Unitigs with no anchor match stay in the raw graph
assembly and are not exported or carried into the next round. Use both
end-anchor flags to preserve the original sequences
through both phases; ordinary flanking rebuilds its graph from whole contigs.
Reduction anchors and raw unitigs are kept under `tmp/kN/assembly_graph/`
when using `--keep-tmp-files`.

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

**Reduction and descending assembly.** For a list such as `21,31,41,51`,
count read edges in ascending order and reduce each smaller k against the next
larger's read-counted edges. All reductions finish before any SdBG is built.
Assemble k51 directly from its counted edges. Then build k41 from the k51 main
contigs and k41 reduced edges, k31 from the k41 main contigs and k31 reduced
edges, and k21 from the k31 main contigs and k21 reduced edges. These steps use
the normal `seq2sdbg` contig input, including contig-header multiplicities, and
standard assembly cleaning. No `merged_edges`, mapping, consensus, or iterative
read-edge extraction is used during the descent.

Reduction uses bounded block I/O, packed-word operations, and chunk-based workers;
its Bloom allocation receives the workflow memory budget. Surviving reduced edges
are exactly deduplicated using buffers sized by `--merge-sort-memory`. See
[Edge Reduction Performance](reduction-performance.md) for parameters and benchmarks.
Assembly skips unused depth inference, uses thread-private histograms, avoids
no-change refreshes, and checks changed-only output before sequence reconstruction.
See [Assembly Performance](assembly-performance.md) for details and comparison commands.
`--merge-index` and `--merge-dedup` are accepted for native CLI and checkpoint
compatibility but no longer affect the workflow. The standalone core merge
remains available. Incomplete runs with older forward-merge checkpoints reuse
eligible count stages and rebuild the changed stages without reading old
`merged_edges`. Completed runs remain completed on `--continue`; use a fresh
output directory to rebuild them with this behavior.

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
