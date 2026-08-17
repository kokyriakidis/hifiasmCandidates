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
    hetmer.cpp
    fakechain.cpp
    fakechain_bridge.cpp
)

# --- Bundled myloasm SNPmer engine (external/myloasm, Rust staticlib). --------
# hetmer.cpp / fakechain*.cpp call its C ABI (myloasm_ffi.h), so the library and
# any consumer must link libmyloasm.a plus the native libs cargo reports. We
# build it with cargo at configure/build time, mirroring the Makefile.
find_program(CARGO_EXECUTABLE cargo)
if(NOT CARGO_EXECUTABLE)
    message(FATAL_ERROR
        "cargo (Rust toolchain) not found. The bundled myloasm SNPmer engine is "
        "required. Install Rust via https://rustup.rs and re-run cmake.")
endif()

set(MYLOASM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/myloasm")
if(NOT EXISTS "${MYLOASM_DIR}/Cargo.toml")
    message(FATAL_ERROR
        "external/myloasm is not checked out. Run:\n"
        "  git submodule update --init --recursive")
endif()
set(MYLOASM_LIB "${MYLOASM_DIR}/target/release/libmyloasm.a")

add_custom_command(
    OUTPUT "${MYLOASM_LIB}"
    COMMAND ${CARGO_EXECUTABLE} build --release
            --manifest-path "${MYLOASM_DIR}/Cargo.toml"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Building bundled myloasm SNPmer engine (cargo build --release)"
    VERBATIM)
add_custom_target(myloasm_engine DEPENDS "${MYLOASM_LIB}")

# Native libs cargo links against for a staticlib (cargo rustc
# --print native-static-libs); required at final link of any consumer.
set(MYLOASM_SYS_LIBS gcc_s util rt dl)

add_library(hifiasm_overlaps STATIC ${HIFIASM_OVERLAPS_SOURCES})
add_dependencies(hifiasm_overlaps myloasm_engine)

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
    # Bundled myloasm staticlib + the native libs it depends on. PUBLIC so a
    # consumer that links hifiasm_overlaps pulls these transitively.
    "${MYLOASM_LIB}"
    ${MYLOASM_SYS_LIBS}
)

# Note: the dead error-correction/assembly code has been pruned from
# Correct.cpp/ecovlp.cpp, so this library links standalone without
# -Wl,--gc-sections. The flag is still applied to our own executable (below)
# to keep the binary minimal, but it is intentionally NOT propagated as an
# INTERFACE requirement so consumers (e.g. dinara's whole-program static
# link) are not forced to enable it.
if(HIFIASM_OVERLAPS_BUILD_EXE)
    add_executable(hifiasm main.cpp)
    target_link_libraries(hifiasm PRIVATE hifiasm_overlaps)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_link_options(hifiasm PRIVATE -Wl,--gc-sections)
    endif()
endif()
