# Random Fuzzy Library - RFL

![Build and Test Status](https://github.com/pauldruce/RFL/actions/workflows/build_and_test_release.yml/badge.svg)

This repository contains the source code for a C++ library to construct and run simulations of Finite Noncommutative Geometries (Finite NCGs).
Note that because this is a library and not an application, the reader must implement their application in C++ and include this project.
Some demo applications are included to demonstrate this.

## Background

Random Finite NCGs are at the forefront of academic research into Finite NCGs, and several interesting results have been found.
Some relevant academic papers can be found here:

1. John W. Barrett and Lisa Glaser. (2016). Monte Carlo simulations of random non-commutative geometries. Journal of Physics A: Mathematical and Theoretical 49, 24: 245001. https://doi.org/10.1088/1751-8113/49/24/245001
2. Lisa Glaser. (2017). Scaling behavior in random non-commutative geometries. Journal of Physics A: Mathematical and Theoretical 50, 27: 275201. https://doi.org/10.1088/1751-8121/aa7424
3. John W Barrett, Paul Druce, and Lisa Glaser. (2019). Spectral estimators for finite non-commutative geometries. Journal of Physics A: Mathematical and Theoretical 52, 27: 275203. https://doi.org/10.1088/1751-8121/ab22f8
4. Lisa Glaser and Abel Stern. (2020). Understanding truncated non-commutative geometries through computer simulations. Journal of Mathematical Physics 61, 3: 033507. https://doi.org/10.1063/1.5131864
5. Lisa Glaser and Abel B. Stern. (2021). Reconstructing manifolds from truncations of spectral triples. Journal of Geometry and Physics. https://doi.org/10.1016/j.geomphys.2020.103921

For an introduction to the area of Non-commutative geometry and specifically finite non-commutative geometry, see:

6. John W. Barrett, "Matrix geometries and fuzzy spaces as finite spectral triples", J. Math. Phys. 56, 082301 (2015) https://doi.org/10.1063/1.4927224

## Dependencies

This project requires you to have the libraries GSL and Armadillo available and CMake installed.

The installation process for both of these projects is very platform dependent, you should follow the installation instructions detailed by each project

- GSL: [https://www.gnu.org/software/gsl/]()
- Armadillo: [http://arma.sourceforge.net]()

On macOS you should be okay using Homebrew, i.e. `brew install gsl` and `brew install armadillo`

On Ubuntu, I imagine `apt-get install libgsl-dev` and `apt-get install libarmadillo-dev` should suffice.

For Windows and other Linux distros you will have to follow the installation instruction on the website above. Good luck :)

## Building the library

Once you have installed the required dependencies, building RFL library can be done in various ways, but the simplest
is to use CMake.
To install CMake, please follow the instructions on the website: [https://cmake.org]()

- On macOS - you can use Homebrew - `brew install cmake`
- On Ubuntu - you can use `apt-get install cmake` should work.
- For Windows and other Linux distros, follow the instruction on the CMake webpage.

Most C/C++ IDEs will have CMake capabilities. However, to build this project manually, you need the following commands

1. Create a directory for your build files within the cloned RFL repo. Typical directory names are `build`, `build-debug`, `build-release` etc.
   A typical set of commands for cloning, and creating the build directory looks like this:
   ```bash
   git clone https://github.com/pauldruce/RFL.git
   cd RFL
   mkdir build
   cd build
   ```
2. Call CMake on the source files. This will create the build files for the project using whatever build system your operating system defaults to.
   Note that it will produce the build files in the current working directory, so make sure to be in a build-only directory before calling Cmake.

   ```bash
   # The command should be "cmake /path/to/source/files" from within an empty build directory
   # In the command below, we are assuming the source files (importantly, the CMakeLists.txt) are
   cmake ..

   # Another handy version of the cmake command is:
   cmake -B ./build .
   ```

3. Build the project. This can be done by manually calling `make` or it's equivalent in the `build` directory. Or CMake has a handy command which is platform independent:
   ```bash
   cmake --build . --target all
   ```
   This build process can be multi-threaded with the CMake command by adding `-j #threads` to the above command. For instance
   ```bash
   cmake --build . --target all -j 4
   ```
   will run the build process on 4 threads, and should be a fair bit quicker. :)
4. Run CTest to check everything works correctly. This is done simply by calling `ctest` from the build directory. Like the `cmake --build` command above,
   the tests can be ran on multiple threads with the `-j` command. For instance, to run the tests on 4 threads, just type:
   ```bash
   ctest -j 4
   ```

## Documentation

Documentation for how this software is architected and implemented can be found in the READMEs in the directory `RFL_source` and its subdirectories.
