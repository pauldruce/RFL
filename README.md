# Random Fuzzy Library - RFL

![Build and Test Status](https://github.com/pauldruce/RFL/actions/workflows/build_and_test.yml/badge.svg)

RFL is a C++ library to construct and run Finite Non-commutative Geometries (Finite NCGs) simulations. See §[Background](#background) below for more details on what a Finite NCG is.

This library is written using C++17 and uses Armadillo and GSL for its mathematical operations.
The library is built using CMake, and detailed instructions on how to generate the library files for your desired platform are listed below in  §[Building the library](#building-the-library)

Note that because this is a library and not an application, the reader must implement their application in C++ and include this project.
Some examples of applications are included in `/examples` to demonstrate how to do this.
These examples can also be built using CMake, and instructions on they are built are in the README in `/examples`

## Library Dependencies
Building this project requires you to have the libraries GSL and Armadillo available and CMake installed.

The installation processes for GSL and Armadillo are very platform dependent, and you should follow the installation instructions detailed by each project:

- GSL: [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)
- Armadillo: [http://arma.sourceforge.net](http://arma.sourceforge.net)

Here are my recommended methods of installing these dependencies per platform:

* On macOS, you should be okay using Homebrew, i.e. `brew install gsl` and `brew install armadillo`
* On Ubuntu, I imagine `apt-get install libgsl-dev` and `apt-get install libarmadillo-dev` should suffice.
* For Windows and other Linux distributions, you must follow the installation instructions on the website above.

To install CMake, please follow the instructions on the website: [https://cmake.org](https://cmake.org)

Here are my recommended methods of installing CMake per platform:

- On macOS - you can use Homebrew - `brew install cmake`
- On Ubuntu - you can use `apt-get install cmake` should work.
- Follow the instructions on the CMake webpage for Windows and other Linux distributions.


## Building the library

Once you have installed the required dependencies, the RFL library can be built using CMake.

tl;dr
```bash
   git clone https://github.com/pauldruce/RFL.git
   cd RFL
   mkdir build
   cd build
   cmake ..
   cmake --build . --target all -j 4
   ctest -j 4
```

Most C/C++ IDEs will have CMake capabilities, and CMake offers a GUI application to make this process easier.

However, to build this project via a terminal, you use the following commands:

1. Create a directory for your build files within the cloned RFL repo. Typical directory names are `build`, `build-debug`, `build-release` etc.
   A typical set of shell commands for cloning this repo and creating the build directory looks like this:
   ```bash
   git clone https://github.com/pauldruce/RFL.git
   cd RFL
   mkdir build
   cd build
   ```
2. __Call CMake on the source files.__ CMake will create the build files for the project using whatever build system your operating system defaults to.
   Note that it will produce the build files in the current working directory, so make sure to be in a build-only directory before calling CMake.
   There are two main ways to do this via the terminal:
   * The command: `cmake ..` from the new `./build` directory created by step 1.
     The command should generally be `cmake /path/to/source/files` from within an empty build directory. We assumed the source files (notably, the CMakeLists.txt) are in the parent directory where we are running this command.
   * The command: `cmake -B ./build .` from the root directory of this project - skipping step 1.
     This is a different but handy way to create the build files. This command will create a folder `./build` and generate the build files within it.
3. __Build the project__ This can be done by manually calling `make`, or it is equivalent in the `build` directory. Or CMake has a handy command which is platform-independent:
   ```bash
   cmake --build . --target all
   ```
   This build process can be multi-threaded with the CMake command by adding `-j #threads` to the above command. For instance,
   ```bash
   cmake --build . --target all -j 4
   ```
   will run the build process on four threads and should be quicker.
4. Run CTest to check everything works correctly. This is done simply by calling `ctest` from the build directory. Similar to the `cmake --build` command above,
   the tests can be run on multiple threads with the `-j` command. For instance, to run the tests on four threads, just type:
   ```bash
   ctest -j 4
   ```



## Documentation

Documentation for how this software is architected and implemented can be found in the READMEs in the
directory `RFL_source` and its subdirectories.


## Background

Random Finite NCGs are at the forefront of academic research into Finite NCGs, and several interesting results have been found.
Some relevant academic papers can be found here:

1. John W. Barrett and Lisa Glaser. (2016). Monte Carlo simulations of random non-commutative geometries. Journal of Physics A: Mathematical and Theoretical 49, 24: 245001. https://doi.org/10.1088/1751-8113/49/24/245001
2. Lisa Glaser. (2017). Scaling behaviour in random non-commutative geometries. Journal of Physics A: Mathematical and Theoretical 50, 27: 275201. https://doi.org/10.1088/1751-8121/aa7424
3. John W Barrett, Paul Druce, and Lisa Glaser. (2019). Spectral estimators for finite non-commutative geometries. Journal of Physics A: Mathematical and Theoretical 52, 27: 275203. https://doi.org/10.1088/1751-8121/ab22f8
4. Lisa Glaser and Abel Stern. (2020). Understanding truncated non-commutative geometries through computer simulations. Journal of Mathematical Physics 61, 3: 033507. https://doi.org/10.1063/1.5131864
5. Lisa Glaser and Abel B. Stern. (2021). Reconstructing manifolds from truncations of spectral triples. Journal of Geometry and Physics. https://doi.org/10.1016/j.geomphys.2020.103921

For an introduction to the area of Non-commutative geometry and specifically finite non-commutative geometry, see:

6. John W. Barrett, "Matrix geometries and fuzzy spaces as finite spectral triples", J. Math. Phys. 56, 082301 (2015) https://doi.org/10.1063/1.4927224