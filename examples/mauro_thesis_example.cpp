//
// Created by Paul Druce on 09/12/2022.
//

#include <iostream>
#include <gsl/gsl_rng.h>
#include <ctime>
#include "Geom24.hpp"

int main() {
  // Initialize the random number generator
  gsl_rng *engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));
  // Clifford module parameters
  int p = 2;
  int q = 0;
  // Matrix algebra dimension
  int n = 32;
  // Coupling constant value
  double g2 = -3;
  // Metropolis scale factor
  double scale = 0.05;
  // Create the Dirac operator
  Geom24 G(p, q, n, g2);
  // Metropolis simulation
  for (int i = 0; i < 100; ++i) {
	// Metropolis evolution for 100 steps
	G.MMC(scale, 100, engine);
	// Print the value of the action
	G.print_S(std::cout);
  }
}