//
// Created by Paul Druce on 09/12/2022.
//

#include "Geom24.hpp"
#include <ctime>
#include <gsl/gsl_rng.h>
#include <iostream>

int main() {
  // Initialize the random number generator
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));
  // Clifford module parameters
  constexpr int p = 2;
  constexpr int q = 0;
  // Matrix algebra dimension
  constexpr int n = 32;
  // Coupling constant value
  constexpr double g2 = -3;
  // Metropolis scale factor
  // Create the Dirac operator
  Geom24 G(p, q, n, g2);
  // Metropolis simulation
  for (int i = 0; i < 100; ++i) {
    constexpr double scale = 0.05;
    // Metropolis evolution for 100 steps
    G.MMC(scale, 100, engine);
    // Print the value of the action
    G.print_S(std::cout);
  }

  gsl_rng_free(engine);
}