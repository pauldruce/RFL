# Example Applications

This directory is the source code for several demo applications to showcase the RFL library.

The source files named `mauro_hmc_example.cpp` and `mauro_thesis_example.cpp` are examples taken from the PhD thesis of Mauro D'Arcangelo - who initially worked on this library ([here](https://github.com/darcangelomauro/RFL)).
These examples use the 'old'/original library implementation Mauro developed and maintained here.

The other example is in the subdirectory `Type13`. This example uses the new/modern implementation of the library.
It demonstrates how to create a simple simulation of random NCGs of Clifford type (1,3).

## Building and Running

These examples will be built if you follow the steps in the main README of this repository.
However, to build and run these examples independently, you can run the following command in the CMake build directory:
`cmake --build . --target mauro_mmc_example mauro_hmc_example Type13Metropolis`.
To run any examples, you must navigate to the executable created. For instance, the mauro_mmc_example executable is at `RFL/{cmake_build_dir}/examples/`.