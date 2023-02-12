//
// Created by Paul Druce on 10/12/2022.
//

#include "Simulation.hpp"
#include "Metropolis.hpp"
#include "DiracOperator.hpp"
#include "Action.hpp"
#include "GslRng.hpp"

int main() {
  double metropolis_scale = 0.2;
  int iter = 10;
  auto rng = std::make_unique<GslRng>();

  DiracOperator dirac(1, 3, 10);
  Action action(-2.7, 1.0);
  auto metropolis = std::make_unique<Metropolis>(metropolis_scale, iter, std::move(rng));
  auto simulation = Simulation(dirac,action,std::move(metropolis));

  for (int i = 0; i < 10; i++) {
    simulation.run();
    action.printS(dirac, std::cout);
  }

  return 0;
}