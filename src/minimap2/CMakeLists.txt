set(MINIMAP2_SOURCES
        main.c kthread.c kalloc.c misc.c bseq.c sketch.c sdust.c options.c
        index.c lchain.c align.c hit.c seed.c jump.c map.c format.c pe.c
        esterr.c splitidx.c ksw2_ll_sse.c)
set(MINIMAP2_SIMD_SOURCES ksw2_extz2_sse.c ksw2_extd2_sse.c ksw2_exts2_sse.c)

add_library(minimap2 STATIC ${MINIMAP2_SOURCES})
set_target_properties(minimap2 PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
target_compile_definitions(minimap2 PRIVATE HAVE_KALLOC)
target_include_directories(minimap2 PRIVATE ${ZLIB_INCLUDE_DIRS})
set_source_files_properties(main.c PROPERTIES COMPILE_DEFINITIONS "main=main_minimap2")

if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i.86)$")
    # Match upstream's runtime-dispatched SSE2/SSE4.1 implementation.
    foreach (variant sse2 sse41)
        add_library(minimap2_${variant} OBJECT ${MINIMAP2_SIMD_SOURCES})
        target_compile_definitions(minimap2_${variant} PRIVATE HAVE_KALLOC KSW_CPU_DISPATCH)
        set_target_properties(minimap2_${variant} PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    endforeach ()
    target_compile_definitions(minimap2_sse2 PRIVATE KSW_SSE2_ONLY)
    target_compile_options(minimap2_sse2 PRIVATE -msse2 -mno-sse4.1)
    target_compile_options(minimap2_sse41 PRIVATE -msse4.1)
    target_compile_options(minimap2 PRIVATE -msse2)
    target_compile_definitions(minimap2 PRIVATE KSW_CPU_DISPATCH)
    target_sources(minimap2 PRIVATE ksw2_dispatch.c
            $<TARGET_OBJECTS:minimap2_sse2> $<TARGET_OBJECTS:minimap2_sse41>)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|arm.*)$")
    target_sources(minimap2 PRIVATE ${MINIMAP2_SIMD_SOURCES})
    target_include_directories(minimap2 PRIVATE sse2neon)
    target_compile_definitions(minimap2 PRIVATE KSW_SSE2_ONLY __SSE2__ _FILE_OFFSET_BITS=64)
    target_compile_options(minimap2 PRIVATE -fsigned-char)
    if (NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$")
        target_compile_options(minimap2 PRIVATE -mfpu=neon)
    endif ()
else ()
    message(FATAL_ERROR "Bundled Minimap2 requires x86 SSE2 or ARM NEON")
endif ()

target_link_libraries(minimap2 PUBLIC ${ZLIB_LIBRARIES} Threads::Threads m)
