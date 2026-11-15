# libsais Integration

Unmodified upstream sources and headers from libsais v2.12.1:

- Repository: https://github.com/IlyaGrebnov/libsais
- Commit: `c284ba843dd16c159acb976f17bc19382793ea75`
- License: Apache License 2.0 with LLVM Exceptions, reproduced in `LICENSE`.
- Included: `libsais.c`, `libsais64.c`, and their public headers.

CMake compiles these C99 sources into a static library with OpenMP support.
The library is linked into the native `megahit` executable and the development
`megahit_core` target. No network download or Python step is needed when building.

`bloom_filter/end_fm_index.cpp` uses the 32-bit suffix-array API when possible and
the 64-bit API for texts beyond the signed 32-bit limit. It constructs a BWT and
rank dictionaries for exact, fixed-length contig-end lookup. This integration
does not change FLASH, Minimap2, read encoding, or the assembly graph format.
