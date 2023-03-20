# RFL Documentation

This directory contains the source code for the RFL library.

There are two main classes defined by the library: `Cliff` and `Geom24` defined in `/include/Cliff.hpp` and `/include/Geom24.hpp` respectively.

`Cliff` is responsible for creating the 'gamma matrices' for a specific Clifford module.
The general way to specify a Clifford module is by specifying two positive integers $p$ and $q$.

`Geom24` is the main Class for this library. It is responsible for setting up and running the simulation.

The code at the moment inherently uses an action of the form:
$$S(D) = g_2* Tr(D^2) + g_4*Tr(D^4) $$
where $g_2$ and $g_4$ are real numbers. Hence the name of the class - `Geom24`
