# C++ API Example

This directory contains a complete C++ application demonstrating how to consume the `RFL::core` library.

The application initialises a Dirac operator with signature $(p=1, q=3)$, runs $100$ Metropolis Monte Carlo update steps with Barrett-Glaser parameters $(g_2=-1.0, g_4=1.0)$, and computes the eigenvalue spectrum.

---

## 1. Build and Run with CMake (Recommended)

From the root directory of the repository:

```bash
# Configure and build all examples
cmake -B build
cmake --build build --target main

# Run the compiled binary
./build/examples/cpp/main
```

---

## 2. Build and Run with Makefile (Linux / Unix / macOS)

A standalone `Makefile` is provided in this directory for direct compilation:

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

You can compile and link `main.cpp` directly with `g++` or `clang++` without using CMake.

### Prerequisites
* Ensure `gsl` and `armadillo` are installed on your system (`brew install gsl armadillo` on macOS, or `sudo apt-get install libgsl-dev libarmadillo-dev` on Debian/Ubuntu).
* Build `librfl_core.a` once using CMake:
  ```bash
  cmake -B build
  cmake --build build --target rfl_core
  ```

### Direct GCC Command
From the root directory of the repository:

```bash
g++ -std=c++17 -O3 examples/cpp/main.cpp \
    -Isrc/core \
    -Lbuild/src/core \
    -lrfl_core -larmadillo -lgsl -lgslcblas \
    -o main_gcc

./main_gcc
```

> [!NOTE]
> If `librfl_core.a` was built with AddressSanitizer enabled (the default in non-Release developer builds), add `-fsanitize=address` to the compile and link flags:
> ```bash
> g++ -std=c++17 -fsanitize=address -O3 examples/cpp/main.cpp \
>     -Isrc/core \
>     -Lbuild/src/core \
>     -lrfl_core -larmadillo -lgsl -lgslcblas \
>     -o main_gcc
> ```

---

## 3. Helper Script

For convenience, run the included shell script:

```bash
./examples/cpp/compile_gcc.sh
```
