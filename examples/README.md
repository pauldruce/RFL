# Example Applications

In this directory is the source code for a number of demo applications, to showcase the RFL library.

The source files named `mauro_hmc_example.cpp` and `mauro_thesis_example.cpp` are examples taken from the PhD thesis of Mauro D'Arcangelo - whose initially worked on this library ([here](https://github.com/darcangelomauro/RFL)). These examples use the 'old'/original library implementation that Mauro developed and that is maintained here.

The other example is housed in the subdirectory `Type13`. This example uses the new/modern implementation of the library.
It demonstrates how to create a simple simulation of random NCGs of Clifford type (1,3).

## Building and Running

These examples will have been built if you followed the steps in the main README of this repository.
However, to build and run these examples on their own, you can run the following command in the cmake build directory:
`cmake --build . --target mauro_mmc_example mauro_hmc_example Type13Metropolis`

To then run any of the examples, you need to navigate to the executable that has been created. For instance, for mauro_mmc_example executable can be found at `RFL/{cmake_build_dir}/examples/`.