//
// Created by Paul Druce on 10/12/2022.
//

#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"
#include "Simulation.hpp"

class Recorder {
public:
  void recordEigenvalues(const DiracOperator& dirac, std::string path = "./output") {
    auto eigenvalues = dirac.getEigenvalues();
    eigenvalues.print();
  }

private:
  std::string m_output_path;
};

int main() {
  double metropolis_scale = 0.2;
  int iter = 10;
  auto rng = std::make_unique<GslRng>();

  auto dirac = std::make_unique<DiracOperator>(1, 3, 10);

  auto action = std::make_unique<Action>(-2.7, 1.0);

  auto metropolis = std::make_unique<Metropolis>(
      std::move(action),
      metropolis_scale,
      iter,
      std::move(rng));

  auto simulation = Simulation(std::move(dirac), std::move(metropolis));

  for (int i = 0; i < 10; i++) {
    simulation.run();
    const auto& dirac2 = simulation.getDiracOperator();
    Recorder recorder;
    recorder.recordEigenvalues(dirac2);
  }

  return 0;
}
