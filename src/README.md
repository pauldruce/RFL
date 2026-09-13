# RFL Source Tree

This directory contains the legacy implementation and the modern C++ core of the Random Fuzzy Library (RFL).

Version 1.0 will deprecate the legacy code. A protected Git branch preserves this historical implementation for baseline comparison.

## Legacy RFL

The legacy source code resides in `src/legacy/`. The library defines two primary classes in `src/legacy/include/`:

* `Cliff` (`Cliff.hpp`): Creates the gamma matrices for a specific Clifford module. Two positive integers, $p$ and $q$, define the Clifford signature.
* `Geom24` (`Geom24.hpp`): Sets up and executes the simulation.

The legacy code evaluates the Barrett-Glaser action:

$$S(D) = g_2 \mathrm{Tr}(D^2) + g_4 \mathrm{Tr}(D^4)$$

where $g_2, g_4 \in \mathbb{R}$. The class name `Geom24` originates from this action, because the action contains quadratic ($D^2$) and quartic ($D^4$) Dirac operator traces.

## Core RFL

The modern library implementation resides in `src/core/`. This modern codebase refactors the legacy code into modular, extensible components without compromising performance.

Key improvements include:
* Value semantics and regular C++ types without hidden pointer ownership.
* Decoupled geometry state, action functionals, and Markov chain steppers.
* Strong type safety and standard C++17 conformance.

To report bugs or propose enhancements, open a GitHub issue.

### Performance Parity Test

The core library preserves the high performance of the legacy implementation. A dedicated benchmark test (`src/core/tests/performance/tBenchmark.cpp`) compares execution times between both implementations. The test requires the modern code to run within 5% of the legacy execution time.
