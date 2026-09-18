# C++ API Example

This directory contains a complete C++ application demonstrating how to consume the `RFL::core` library.

The application initialises a Dirac operator with signature $(p=1, q=3)$. It executes $100$ Metropolis update steps with parameters $(g_2=-1.0, g_4=1.0)$ and computes the eigenvalue spectrum.

---

## 1. Build and Run with CMake (Recommended)

From the repository root directory:

```bash
# Configure and build all examples
cmake -B build
cmake --build build --target main

# Run the compiled binary
./build/examples/cpp/main
```

---

## 2. Build and Run with Makefile (Linux / Unix / macOS)

This directory provides a standalone `Makefile` for direct compilation:

```bash
cd examples/cpp

# Build the example binary
make

# Run the simulation
make run

# Clean build artifacts
make clean
```

---

## 3. Compile and Link Directly with GCC or Clang

Compile and link `main.cpp` directly with `g++` or `clang++` without CMake.

### Prerequisites
* Install `gsl` and `armadillo` (`brew install gsl armadillo` on macOS, or `sudo apt-get install libgsl-dev libarmadillo-dev` on Debian/Ubuntu).
* Build `librfl_core.a` once using CMake:
  ```bash
  cmake -B build
  cmake --build build --target rfl_core
  ```

### Direct GCC Command
From the repository root directory:

```bash
mkdir -p build/examples/cpp
g++ -std=c++17 -O3 examples/cpp/main.cpp \
    -Isrc/core \
    -Lbuild/src/core \
    -lrfl_core -larmadillo -lgsl -lgslcblas \
    -o build/examples/cpp/main_gcc

./build/examples/cpp/main_gcc
```

> [!NOTE]
> When building `librfl_core.a` with AddressSanitizer enabled (`-DRFL_ENABLE_ASAN=ON`), pass `-fsanitize=address` to both compiler and linker flags:
> ```bash
> mkdir -p build/examples/cpp
> g++ -std=c++17 -fsanitize=address -O3 examples/cpp/main.cpp \
>     -Isrc/core \
>     -Lbuild/src/core \
>     -lrfl_core -larmadillo -lgsl -lgslcblas \
>     -o build/examples/cpp/main_gcc
> ```

---

## 4. Helper Script

For convenience, execute the included shell script:

```bash
./examples/cpp/compile_gcc.sh
```
