#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "DiracOperator.hpp"
#include "StdRng.hpp"
#include <iostream>
#include <memory>

int main() {
  auto dirac = std::make_unique<DiracOperator>(1, 3, 6);
  auto action = std::make_unique<Action>(-1.0, 1.0);
  auto rng = std::make_unique<StdRng>(42UL);

  Metropolis metropolis(std::move(action), 0.1, 10, std::move(rng));
  double acceptance = metropolis.updateDirac(*dirac);

  std::cout << "Downstream Canary Dirac Type: (" << dirac->getType().first << ", "
            << dirac->getType().second << ")\n";
  std::cout << "Downstream Canary Acceptance Rate: " << acceptance * 100.0 << "%\n";

  return 0;
}
