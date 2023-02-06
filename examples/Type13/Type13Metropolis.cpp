//
// Created by Paul Druce on 10/12/2022.
//

#include "Simulation.hpp"
#include "Metropolis.hpp"
#include "DiracOperator.hpp"
#include "Action.hpp"

int main() {
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  double metropolis_scale = 0.2;
  int iter = 10;

  DiracOperator dirac(1, 3, 10);
  Action action(-2.7, 1.0);
  auto metropolis = Metropolis(metropolis_scale, iter, engine);
  auto simulation = Simulation(dirac,action,metropolis);

  for (int i = 0; i < 10; i++) {
    simulation.run();
    action.printS(dirac, std::cout);
  }

  gsl_rng_free(engine);
  return 0;
}