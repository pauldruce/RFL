//
// Created by Paul Druce on 10/12/2022.
//

#include "Simulation.hpp"
#include "Metropolis.hpp"
#include "DiracOperator.hpp"
#include "Action.hpp"

int main() {
  DiracOperator D(1, 3, 10);
  Action A(-2.7, 1.0);
  Metropolis M;
  Simulation Sim(D, A, M);

  const int stepSize = 100;
  const double scale = 0.05;
  for (int i = 0; i < 10; i++) {
	Sim.run(scale, stepSize);
	A.print_S(D, std::cout);
  }
  return 0;
}