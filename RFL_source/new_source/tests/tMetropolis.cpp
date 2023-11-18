//
// Created by Paul Druce on 04/02/2023.
//

#include "BarrettGlaser/Metropolis.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"
#include <gtest/gtest.h>

TEST(MetropolisTests, ConstructorDoesNotThrow) {
  // TODO: How do we test random stuff? -> find out
  auto rng = std::make_unique<GslRng>();
  constexpr double scale = 0.1;
  constexpr int num_steps = 20;
  auto action = std::make_unique<Action>(1.0, 1.0);

  ASSERT_NO_THROW(
      Metropolis metropolis(std::move(action), scale, num_steps, std::move(rng)););
}

TEST(MetropolisTests, UpdateDiracUpdatesTheDirac) {
  auto rng = std::make_unique<GslRng>();
  constexpr double scale = 0.2;
  constexpr int num_steps = 20;
  auto action = std::make_unique<Action>(1.0, 1.0);

  const Metropolis metropolis(std::move(action), scale, num_steps, std::move(rng));

  const auto dirac = DiracOperator(1, 1, 5);
  const auto old_dirac_matrix = dirac.getDiracMatrix();
  metropolis.updateDirac(dirac);
  const auto new_dirac_matrix = dirac.getDiracMatrix();

  const auto diracs_are_equal = arma::approx_equal(new_dirac_matrix, old_dirac_matrix, "absdiff", 1e-6);

  ASSERT_FALSE(diracs_are_equal) << "The dirac matrix should have been changed";
}