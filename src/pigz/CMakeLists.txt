cmake_minimum_required(VERSION 3.5...4.1)
project(megahit_pigz C)

option(STATIC_BUILD "Link bundled pigz statically" OFF)
if(STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

add_executable(pigz pigz.c yarn.c try.c)
set_target_properties(pigz PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}")
target_compile_definitions(pigz PRIVATE NOZOPFLI _FILE_OFFSET_BITS=64)
target_include_directories(pigz PRIVATE ${ZLIB_INCLUDE_DIRS})
target_link_libraries(pigz PRIVATE ${ZLIB_LIBRARIES} Threads::Threads m)
if(STATIC_BUILD)
    set_target_properties(pigz PROPERTIES LINK_FLAGS "-static")
endif()
