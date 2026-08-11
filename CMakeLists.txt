# CMake build for the stripped-down hifiasm candidate-overlap detector.
#
# Targets:
#   hifiasm_overlaps  - static library exposing the C API in hifiasm_overlaps.h
#                       (built PIC so it can also be linked into shared objects).
#   hifiasm           - CLI executable (optional, ON by default).
#
# The library is meant to be consumed as a submodule, e.g. from dinara:
#     add_subdirectory(external/hifiasmCandidates)
#     target_link_libraries(dinara PRIVATE hifiasm_overlaps)
# The imported target already carries the include directory, the required
# system libraries (z, pthread, m) and the --gc-sections link flag needed to
# drop the unused error-correction/assembly code that the overlap path
# references but never calls.

cmake_minimum_required(VERSION 3.13)
project(hifiasm_overlaps LANGUAGES C CXX)

option(HIFIASM_OVERLAPS_BUILD_EXE "Build the hifiasm CLI executable" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)

# Position-independent code so the static lib can be linked into executables or
# shared libraries alike.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# All translation units except main.cpp (the CLI entry point).
set(HIFIASM_OVERLAPS_SOURCES
    CommandLines.cpp
    Process_Read.cpp
    Hash_Table.cpp
    kthread.cpp
    htab.cpp
    hist.cpp
    sketch.cpp
    anchor.cpp
    extract.cpp
    sys.cpp
    kalloc.cpp
    Correct.cpp
    ecovlp.cpp
    candidates.cpp
    hifiasm_overlaps.cpp
)

add_library(hifiasm_overlaps STATIC ${HIFIASM_OVERLAPS_SOURCES})

# Match the Makefile's architecture/optimization flags. These ISA flags target
# x86-64; adjust for other architectures.
target_compile_options(hifiasm_overlaps PRIVATE
    -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall
    -ffunction-sections -fdata-sections
)

# Consumers get the public header directory automatically.
target_include_directories(hifiasm_overlaps PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(hifiasm_overlaps PUBLIC
    ZLIB::ZLIB Threads::Threads m
)

# The overlap code path references (but never calls) error-correction/assembly
# symbols in Correct.cpp/ecovlp.cpp. --gc-sections lets the linker drop those
# unused sections; propagate it to anything that links this library so their
# final link succeeds without the missing symbols.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_options(hifiasm_overlaps INTERFACE -Wl,--gc-sections)
endif()

if(HIFIASM_OVERLAPS_BUILD_EXE)
    add_executable(hifiasm main.cpp)
    target_link_libraries(hifiasm PRIVATE hifiasm_overlaps)
endif()
