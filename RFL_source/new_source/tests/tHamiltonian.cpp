//
// Created by Paul Druce on 04/02/2023.
//

#include "BarrettGlaser/Hamiltonian.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"
#include <gtest/gtest.h>

TEST(HamiltonianTests, ConstructorDoesNotThrow) {
  constexpr Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  constexpr double step_size = 0.1;
  auto action = std::make_unique<Action>(1.0, 1.0);
  ASSERT_NO_THROW(
      const Hamiltonian hamiltonian(std::move(action), integrator, step_size, std::make_unique<GslRng>()););
}

TEST(HamiltonianTests, CanChangeIntegrator) {
  constexpr Integrator integrator = LEAPFROG;
  constexpr double step_size = 0.1;
  auto action = std::make_unique<Action>(1.0, 1.0);
  Hamiltonian hamiltonian(std::move(action), integrator, step_size, std::make_unique<GslRng>());

  ASSERT_EQ(hamiltonian.getIntegrator(), integrator);

  constexpr Integrator new_integrator = OMELYAN;
  hamiltonian.setIntegrator(new_integrator);
  ASSERT_EQ(hamiltonian.getIntegrator(), new_integrator);
}

TEST(HamiltonianTests, CanChangeStepSize) {
  constexpr Integrator integrator = LEAPFROG;
  constexpr double step_size = 0.1;
  auto action = std::make_unique<Action>(1.0, 1.0);
  Hamiltonian hamiltonian(std::move(action), integrator, step_size, std::make_unique<GslRng>());

  ASSERT_EQ(hamiltonian.getStepSize(), step_size);

  constexpr double new_step_size = 0.5;
  hamiltonian.setStepSize(new_step_size);
  ASSERT_EQ(hamiltonian.getStepSize(), new_step_size);
}

TEST(HamiltonianTests, UpdateDiracUpdatesTheDirac) {
  const Hamiltonian hamiltonian(std::make_unique<Action>(), Integrator::LEAPFROG, 0.2, std::make_unique<GslRng>());
  const auto dirac = DiracOperator(1, 1, 5);
  const auto old_dirac_matrix = dirac.getDiracMatrix();
  hamiltonian.updateDirac(dirac);

  const auto new_dirac_matrix = dirac.getDiracMatrix();

  const auto diracs_are_equal = arma::approx_equal(new_dirac_matrix, old_dirac_matrix, "absdiff", 1e-6);

  ASSERT_FALSE(diracs_are_equal) << "The dirac matrix should have been changed";
}