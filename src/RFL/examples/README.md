# Example Applications

This directory contains example applications for the RFL library in C++ and Python.

## Directory Structure

- **`cpp/`**: Basic C++ API example. Initialises a Dirac operator, runs the Metropolis algorithm, and computes eigenvalues.
- **`python/`**: Python API example. Demonstrates Python bindings, Metropolis sampling, and NumPy integration.
- **`case_studies/`**: Specialised simulations, parameter tuning, and historical thesis models.

## Build and Execution

Build the examples using CMake in your build directory:

```bash
cmake --build .
```

To run the basic C++ example:

```bash
./src/RFL/examples/cpp/main
```

To run the Python example:

```bash
python3 src/RFL/examples/python/main.py
```