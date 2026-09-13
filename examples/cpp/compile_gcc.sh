#!/usr/bin/env bash
# ==============================================================================
# Compile and run examples/cpp/main.cpp directly with GCC/Clang
# ==============================================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

CXX="${CXX:-g++}"

# Ensure librfl_core.a is built
if [ ! -f "build/src/core/librfl_core.a" ]; then
    echo "Building librfl_core.a first via CMake..."
    cmake -B build
    cmake --build build --target rfl_core
fi

# Detect whether librfl_core was built with AddressSanitizer
SAN_FLAGS=""
if nm build/src/core/librfl_core.a 2>/dev/null | grep -q "___asan"; then
    SAN_FLAGS="-fsanitize=address"
fi

echo "Compiling examples/cpp/main.cpp using ${CXX}..."
${CXX} -std=c++17 ${SAN_FLAGS} -O3 examples/cpp/main.cpp \
    -Isrc/core \
    -Lbuild/src/core \
    -lrfl_core -larmadillo -lgsl -lgslcblas \
    -o examples/cpp/main_gcc

echo "Running examples/cpp/main_gcc..."
./examples/cpp/main_gcc
