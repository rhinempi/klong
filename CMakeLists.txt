cmake_minimum_required(VERSION 3.5...4.1)
project(megahit)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
    set(CMAKE_CXX_STANDARD 11)
endif ()

option(COVERAGE "Generate coverage report" OFF)
option(STATIC_BUILD "Build static executation" OFF)
option(SANITIZER "Enable sanitizers" OFF)
option(TSAN "Enable thread sanitizers" OFF)

include_directories(src)
FILE(GLOB_RECURSE ASMBL_SOURCE "src/assembly/*.cpp")
FILE(GLOB_RECURSE LCASM_SOURCE "src/localasm/*.cpp")
FILE(GLOB_RECURSE IDBA_SOURCE "src/idba/*.cpp")
FILE(GLOB_RECURSE SDBG_SOURCE "src/sdbg/*.cpp")
FILE(GLOB_RECURSE CX1_SOURCE "src/sorting/*.cpp")
FILE(GLOB_RECURSE SEQ_SOURCE "src/sequence/*.cpp")
FILE(GLOB_RECURSE TOOLKIT_SOURCE "src/tools/*.cpp")
FILE(GLOB_RECURSE BLOOM_FILTER_SOURCE "bloom_filter/*.cpp")

LIST(APPEND OTHER_SOURCE
        src/main.cpp
        src/main_assemble.cpp
        src/main_buildlib.cpp
        src/main_iterate.cpp
        src/main_local_assemble.cpp
        src/main_sdbg_build.cpp
        src/iterate/contig_end_consensus.cpp
        src/utils/options_description.cpp
        )

if (STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif (STATIC_BUILD)

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

# Keep the upstream C implementation separate from MEGAHIT's C++ sources.
add_library(flash STATIC
        src/flash/flash.c
        src/flash/combine_reads.c
        src/flash/iostream.c
        src/flash/read_io.c
        src/flash/read_queue.c
        src/flash/read_util.c
        src/flash/util.c)
set_target_properties(flash PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
target_compile_definitions(flash PRIVATE _GNU_SOURCE _FILE_OFFSET_BITS=64)
target_include_directories(flash PRIVATE ${ZLIB_INCLUDE_DIRS})
set_source_files_properties(src/flash/flash.c PROPERTIES
        COMPILE_DEFINITIONS "main=main_flash")
target_link_libraries(flash PUBLIC ${ZLIB_LIBRARIES} Threads::Threads)
add_subdirectory(src/minimap2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXXH_INLINE_ALL -ftemplate-depth=3000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprefetch-loop-arrays -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XROOT__='\"${CMAKE_SOURCE_DIR}/src\"'")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__XFILE__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ZLIB_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O1 -D_LIBCPP_DEBUG -D_GLIBCXX_DEBUG")

if (COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
    set(COV_PY "coverage run")
endif (COVERAGE)

if (SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
endif (SANITIZER)

if (TSAN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
endif (TSAN)


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS}")

set(MEGAHIT_CORE_SOURCE
        ${OTHER_SOURCE}
        ${ASMBL_SOURCE}
        ${IDBA_SOURCE}
        ${SDBG_SOURCE}
        ${LCASM_SOURCE}
        ${SEQ_SOURCE}
        ${CX1_SOURCE}
        ${TOOLKIT_SOURCE}
        ${BLOOM_FILTER_SOURCE})

add_executable(megahit_core ${MEGAHIT_CORE_SOURCE})


if (STATIC_BUILD)
    # TODO dirty
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -static")
    set_target_properties(megahit_core PROPERTIES LINK_SEARCH_START_STATIC ON)
    set_target_properties(megahit_core PROPERTIES LINK_SEARCH_END_STATIC ON)
endif (STATIC_BUILD)

target_link_libraries(megahit_core ${ZLIB_LIBRARIES} flash minimap2)

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/megahit
        COMMAND python3 ${CMAKE_SOURCE_DIR}/cmake/embed_megahit.py
                --launcher ${CMAKE_SOURCE_DIR}/src/megahit
                --core $<TARGET_FILE:megahit_core>
                --output ${CMAKE_CURRENT_BINARY_DIR}/megahit
        DEPENDS megahit_core ${CMAKE_SOURCE_DIR}/src/megahit ${CMAKE_SOURCE_DIR}/cmake/embed_megahit.py
        VERBATIM)

add_custom_target(embedded_megahit ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/megahit)

set(TEST_DATA ${CMAKE_SOURCE_DIR}/test_data)

add_custom_target(
        simple_test
        COMMAND ./megahit --test -t 2 -o test-embedded --force
)

add_dependencies(simple_test embedded_megahit)

add_custom_target(flash_test
        COMMAND ${CMAKE_COMMAND} -E env
                "MEGAHIT_TEST_CORE=$<TARGET_FILE:megahit_core>"
                python3 -m unittest discover -s ${CMAKE_SOURCE_DIR}/tests -p test_flash.py
        DEPENDS megahit_core
        VERBATIM)

add_custom_target(minimap2_test
        COMMAND ${CMAKE_COMMAND} -E env
                "MEGAHIT_TEST_CORE=$<TARGET_FILE:megahit_core>"
                python3 -m unittest discover -s ${CMAKE_SOURCE_DIR}/tests -p test_minimap2.py
        DEPENDS megahit_core
        VERBATIM)

add_custom_target(contig_ends_test
        COMMAND python3 -m unittest discover -s ${CMAKE_SOURCE_DIR}/tests -p test_contig_ends.py
        VERBATIM)

add_custom_target(end_consensus_test
        COMMAND ${CMAKE_COMMAND} -E env
                "MEGAHIT_TEST_CORE=$<TARGET_FILE:megahit_core>"
                python3 -m unittest discover -s ${CMAKE_SOURCE_DIR}/tests -p test_end_consensus.py
        DEPENDS megahit_core
        VERBATIM)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/megahit DESTINATION bin)
install(DIRECTORY test_data DESTINATION share/${PROJECT_NAME})
