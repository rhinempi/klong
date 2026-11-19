MEGAHIT
=======

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/megahit.svg?style=flat-square&label=BioConda%20install)](https://anaconda.org/bioconda/megahit)
[![Downloads](https://img.shields.io/github/downloads/voutcn/megahit/total?style=flat-square)](https://github.com/voutcn/megahit/releases)
[![Build Status](https://img.shields.io/travis/voutcn/megahit?style=flat-square)](https://travis-ci.org/voutcn/megahit)
[![codecov](https://img.shields.io/codecov/c/github/voutcn/megahit?style=flat-square)](https://codecov.io/gh/voutcn/megahit)

MEGAHIT is an ultra-fast and memory-efficient NGS assembler. It is optimized for metagenomes, but also works well on generic single genome assembly (small or mammalian size) and single-cell assembly.

This fork builds a **native single-file `megahit` executable** with a C++ workflow
driver and linked assembly core. Python is not required to build or run it.
See [native executable details](docs/native-executable.md) for dependencies,
compression, checkpoints, and migration from the historical Python driver.

Installation
---------------

### Conda
```sh
conda install -c bioconda megahit
```

### Guix
```sh
guix install megahit
```

### Pre-built binaries for x86_64 Linux

```sh
wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
cd MEGAHIT-1.2.9-Linux-x86_64-static/bin/
./megahit --test  # run on a toy dataset
./megahit -1 MY_PE_READ_1.fq.gz -2 MY_PE_READ_2.fq.gz -o MY_OUTPUT_DIR
```

### Pre-built docker image
``` sh
# in the directory with the input reads
docker run -v $(pwd):/workspace -w /workspace --user $(id -u):$(id -g) vout/megahit \
  -1 MY_PE_READ_1.fq.gz -2 MY_PE_READ_2.fq.gz -o MY_OUTPUT_DIR
```

### Building from source

#### Prerequisites

-   For building: zlib and bzip2 development headers, CMake &gt;= 3.5, and C/C++11 compilers with OpenMP.
-   For running: standard OS shared libraries (C/C++ runtime, OpenMP, zlib, bzip2). FLASH, Minimap2, and optional pigz are bundled. Python is only needed for the optional regression-test harness.

On Ubuntu/Debian, install the prerequisites with:

```sh
sudo apt-get update
sudo apt-get install -y git build-essential cmake zlib1g-dev libbz2-dev
```

```sh
git clone --branch master https://github.com/rhinempi/klong.git
cd klong
git log -1 --oneline
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -- -j4
./megahit --version
# Optional toy-data validation:
cmake --build . --target simple_test
# Optional installation to a user-owned directory:
# Add -DCMAKE_INSTALL_PREFIX="$HOME/.local" to the configure command, then:
# cmake --build . --target install
```

The executable is `build/megahit`. Pigz is compiled and embedded by default;
if its optional build fails, the executable still builds and uses native zlib gzip.
The build commands above and pigz helper also work with CMake older than 3.13.
No submodule initialization or separate FLASH, Minimap2, or pigz installation is needed.

Usage
-----

### Basic usage
```sh
megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
megahit contig2fastg 119 out/intermediate_contigs/flanking/k119.contigs.fa > k119.fastg
```
The contigs can be found `final.contigs.fa` in the output directory.

### Advanced usage
- Merged FASTA compression uses bundled C pigz with `-t` compression threads and gzip level 1. Pigz source is vendored in `src/pigz` and built/embedded automatically; no system installation is needed. If its optional build fails, or `-DBUNDLED_PIGZ=OFF` is set, MEGAHIT still builds and uses native C zlib gzip. `-t 1` also uses zlib. System PATH is not searched. FASTA formatting streams without an uncompressed intermediate file; backend, thread count, and conversion time are logged. See [pigz build details](src/pigz/INTEGRATION.md).
- Native reduction assemblies retain only contigs longer than the next k (k21 -> k29 retains at least 30 bp). Merging extracts packed (next-k+1)-mer windows directly from these contigs and combines them with newly counted next-k edges; short-contig overlap searching and its hash/FM index are skipped. Intermediate merges defer deduplication until after reduction: only surviving edges are sorted and deduplicated into `reduced_edges`, retaining maximum multiplicity. The last reduction k still deduplicates during merging because no later reduction exists. Extraction, buffer sorting, and key-range merging are parallelized. The final reduction assembly consumes the last merged edge set and keeps its usual output lengths. See [native workflow details](docs/native-executable.md).
- String k-mer/edge reports are disabled in the launcher, including resumed runs. Binary `.edges.*` shards and required `.edges.info` metadata are retained. Reduction does not allocate its diagnostic origin hash map, decode longer edges for provenance, or collect report buffers/counters. Native debugging remains opt-in with `prefixbf --human_readable 1` or `MEGAHIT_HUMAN_READABLE_OUTPUT=1` for core edge writers.
- Paired FASTQ inputs run through bundled [FLASH](https://github.com/ebiggers/flash) before binary-library construction by default. `reads.lib.bin` uses merged reads plus both unmerged mate outputs as single-end sequences; successfully merged pairs are not counted again as original mates. Consensus mapping uses the same merged/unmerged exports. Files are kept in `OUT_DIR/tmp/flash/`, with nonempty inputs listed in `assembly_inputs.json`. Supports `-1`/`-2`, `--12`, gzip, and bzip2 inputs. Use `--no-merge-pairs` for original mate-based input, including paired FASTA. See [integration details](src/flash/INTEGRATION.md).
- After reduction/merge assembly through k59 (default list: 21,29,39,59), consensus extensions seed a fresh k21 graph with retained k59 contigs, bubbles, and additional contigs. The flanking list defaults to 21,31,41,61,81,101,141; both phase lists are independently configurable. Contigs are stored in `intermediate_contigs/flanking/`, and graphs, contig ends, and mappings in `tmp/flanking/kN/`, preserving first-pass outputs. Each transition maps merged and unmerged FLASH exports with Minimap2 (`-x sr -a -Y --secondary=no -t NUM_CPU_THREADS -I 55G`) and publishes only extended consensus records as `kN.local.fa` plus `.info`, with read-support-based multiplicity. `--no-merge-pairs` disables FLASH-based consensus; `--no-local` disables local contributions. See [consensus details](docs/contig-end-consensus.md).
- `tmp/k59/k59.contig_ends.fa` contains the first and last 290 bases of each current contig. Names end in `|start` or `|end`, with 1-based inclusive source coordinates. Contigs shorter than 580 bases retain their complete sequence and original header; exactly 580 bases are split. Ending fragments keep their original orientation. These intermediate end files are not graph inputs themselves; only extended consensus records are published as local-contig inputs. End FASTA, extended FASTA/info, and mapping results survive normal temporary-file cleanup.
- `--kmin-1pass`: if sequencing depth is low and too much memory used when build the graph of k_min
- `--presets meta-large`: if the metagenome is complex (i.e., bio-diversity is high, for example soil metagenomes)
- `--cleaning-rounds 1 --disconnect-ratio 0`: get less pruned assembly (usually shorter contigs)
- `--continue -o out`: resume an interrupted job from `out`

To see the full manual, run `megahit` without parameters or with `-h`.

Also, our [wiki](https://github.com/voutcn/megahit/wiki) may be helpful.

### Examples

There is an example of a real assembly on this [wiki page](https://github.com/voutcn/megahit/wiki/An-example-of-real-assembly) 
and its corresponding [interactive Jupyter notebook](https://biotutorials.org/megahit). Contributed by Noah Lebovic (@lebovic).

Publications
------------

-   Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, doi: 10.1093/bioinformatics/btv033 \[PMID: [25609793](http://www.ncbi.nlm.nih.gov/pubmed/25609793)\].
-   Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.

License
-------

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
