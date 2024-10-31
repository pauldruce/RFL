//
// Created by Paul Druce on 10/12/2022.
//

#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "DiracOperator.hpp"
#include "EigenvalueRecorder.hpp"
#include "GslRng.hpp"
#include "Simulation.hpp"
#include <iomanip>

using namespace arma;

static std::string getCurrentDateTime() {
  const auto t = std::time(nullptr);
  const auto tm = *std::localtime(&t);
  std::stringstream ss;
  ss << std::put_time(&tm, "%Y%m%d%H%M%S");
  return ss.str();
}

int main() {
  double metropolisScale = 0.2;
  int iter = 10;
  auto rng = std::make_unique<GslRng>();
  auto dirac = std::make_unique<DiracOperator>(1, 3, 10);

  auto g2 = -2.7;
  auto g4 = 1.0;
  auto action = std::make_unique<Action>(g2, g4);

  auto metropolis = std::make_unique<Metropolis>(
      std::move(action),
      metropolisScale,
      iter,
      std::move(rng));

  const auto simulation = Simulation(std::move(dirac), std::move(metropolis));
  const auto timeStamp = getCurrentDateTime();

  for (int i = 0; i < 10; i++) {
    simulation.run();
    const auto& dirac2 = simulation.getDiracOperator();
    EigenvalueRecorder recorder(dirac2, g2, timeStamp);
    recorder.recordEigenvalues(i);
  }

  return 0;
}
