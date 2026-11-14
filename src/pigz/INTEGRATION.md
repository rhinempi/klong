# Bundled pigz

Upstream: https://www.zlib.net/pigz/
Version: 2.8 (19 August 2023).
Source archive: https://www.zlib.net/pigz/pigz-2.8.tar.gz
SHA256: `eb872b4f0e1f0ebe59c9f7bd8c506c4204893ba6a8492de31df416f0d5170fd0`.

The upstream C sources and license notices are unchanged. Licensing is the
permissive zlib-style license in the headers of pigz.c, yarn.c, and try.c.
Only the standard compression implementation is built (`NOZOPFLI`); optional
level-11 Zopfli sources are not required for the workflow's level-1 gzip.

With `BUNDLED_PIGZ=ON` (default), packaging builds pigz in an isolated CMake
sub-build and embeds the executable alongside the MEGAHIT core. No download,
system installation, administrator privileges, or system pigz are required.
The sub-build uses the main C compiler and zlib paths. `STATIC_BUILD=ON` is
also forwarded. Its source is vendored so normal builds work offline.

Compiler, linker, or configuration failure emits a warning but does not fail
MEGAHIT packaging. A stale executable is discarded before trying the build.
The resulting package then uses Python gzip, even if PATH contains pigz.
Use `-DBUNDLED_PIGZ=OFF` to deliberately omit the optional compressor.

The self-contained executable extracts the bundled pigz next to its core in
the runtime cache. For `-t > 1`, FASTA compression selects only that bundled
executable. For `-t 1`, or when the bundle is absent/unavailable, Python gzip
is used. System PATH is not searched. A runtime compression error remains an
error, not a silent fallback; partial output is removed and not published.
