# RFL Documentation

This directory contains the source code for the original RFL library created by Mauro D'arcangelo and its modern
implementation.

The aim is to deprecate Mauro's implementation (Original RFL) for version 1.0 of this project.
The code will be preserved under a protected branch on this repository.

## Original RFL

The source code for the original library is in `./original_source`.
The library defines only two classes: `Cliff` and `Geom24`, defined in `/include/Cliff.hpp` and `/include/Geom24.hpp`,
respectively.

- `Cliff` is responsible for creating the 'gamma matrices' for a specific Clifford module.
  The general way to specify a Clifford module is by setting two positive integers $p$ and $q$.
- `Geom24` is the main Class for this library. It is responsible for setting up and running the simulation.

The original RFL code inherently uses an action of the form:
$$S(D) = g_2* Tr(D^2) + g_4*Tr(D^4) $$
where $g_2$ and $g_4$ are real numbers. This action is the origin of the name of the class `Geom24` as the action
contains the quadratic (D^2) and quartic (D^4) traces of the Dirac operator.

## New RFL

The new implementation of RFL has taken the original RFL source code and made it more flexible without compromising
performance. A significant part of the refactoring was to break down the two behemoth classes into smaller classes,
implement a shallow layer of inheritance of abstract classes and migrate memory management to modern C++ techniques (
read smart pointers with C++17).
The new source code is found in `./new_source`.

This refactoring allows for better testing and more flexibility for library users.
It is constantly being developed, and any improvements or enhancements are welcome. Please leave an issue on the GitHub
repo.

The high performance that Mauro's original implementation strived to develop has been preserved.
To ensure that this performance does not degrade, a unit test that captures the difference between the original and new
implementations.
This test requires that the new implementation is no more than 5% slower to run the same simulation in each
implementation.
