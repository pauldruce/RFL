//
// Created by Paul Druce on 04/02/2023.
//

#include "Hamiltonian.hpp"
#include <gtest/gtest.h>

TEST(HamiltonianTests, ConstructorDoesNotThrow) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  double step_size = 0.1;
  ASSERT_NO_THROW(
      const Hamiltonian H(integrator, engine, step_size););
}

TEST(HamiltonianTests, CanChangeEngine) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  double step_size = 0.1;
  Hamiltonian H(integrator, engine, step_size);

  ASSERT_EQ(H.getEngine(), engine);

  const gsl_rng* new_engine = gsl_rng_alloc(gsl_rng_default);
  H.setEngine(new_engine);
  ASSERT_EQ(H.getEngine(), new_engine);
}

TEST(HamiltonianTests, CanChangeIntegrator) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  double step_size = 0.1;
  Hamiltonian H(integrator, engine, step_size);

  ASSERT_EQ(H.getIntegrator(), integrator);

  const Integrator new_integrator = OMELYAN;
  H.setIntegrator(new_integrator);
  ASSERT_EQ(H.getIntegrator(), new_integrator);
}

TEST(HamiltonianTests, CanChangeStepSize) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  double step_size = 0.1;
  Hamiltonian H(integrator, engine, step_size);

  ASSERT_EQ(H.getStepSize(), step_size);

  const double new_step_size = 0.5;
  H.setStepSize(new_step_size);
  ASSERT_EQ(H.getStepSize(), new_step_size);
}

TEST(HamiltonianTests, UpdateDiracUpdatesTheDirac) {
  Hamiltonian H(Integrator::LEAPFROG, gsl_rng_alloc(gsl_rng_ranlxd1), 0.2);
  auto D = DiracOperator(1, 1, 5);
  auto old_D_mat = D.getDiracMatrix();
  auto A = Action(-2.7);
  H.updateDirac(D, A);

  auto new_D_mat = D.getDiracMatrix();

  const auto diracsAreEqual = arma::approx_equal(new_D_mat, old_D_mat, "absdiff", 1e-6);

  ASSERT_FALSE(diracsAreEqual) << "The dirac matrix should have been changed";
}