//
// Created by Paul Druce on 04/02/2023.
//

#include "Hamiltonian.hpp"
#include <gtest/gtest.h>

class Engine {
public:
  Engine(){
    rng = gsl_rng_alloc(gsl_rng_ranlxd1);
  }
  ~Engine(){
    if(rng) {
      gsl_rng_free(rng);
    }
  }
  gsl_rng* rng;
};

TEST(HamiltonianTests, ConstructorDoesNotThrow) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const Engine engine;
  double step_size = 0.1;
  ASSERT_NO_THROW(
      const Hamiltonian hamiltonian(integrator, engine.rng, step_size););
}

TEST(HamiltonianTests, CanChangeEngine) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const Engine engine;
  double step_size = 0.1;
  Hamiltonian hamiltonian(integrator, engine.rng, step_size);

  ASSERT_EQ(hamiltonian.getEngine(), engine.rng);

  const Engine new_engine;
  hamiltonian.setEngine(new_engine.rng);
  ASSERT_EQ(hamiltonian.getEngine(), new_engine.rng);
}

TEST(HamiltonianTests, CanChangeIntegrator) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const Engine engine;

  double step_size = 0.1;
  Hamiltonian hamiltonian(integrator, engine.rng, step_size);

  ASSERT_EQ(hamiltonian.getIntegrator(), integrator);

  const Integrator new_integrator = OMELYAN;
  hamiltonian.setIntegrator(new_integrator);
  ASSERT_EQ(hamiltonian.getIntegrator(), new_integrator);
}

TEST(HamiltonianTests, CanChangeStepSize) {
  const Integrator integrator = LEAPFROG;
  // TODO: How do we test random stuff? -> find out
  const Engine engine;
  double step_size = 0.1;
  Hamiltonian hamiltonian(integrator, engine.rng, step_size);

  ASSERT_EQ(hamiltonian.getStepSize(), step_size);

  const double new_step_size = 0.5;
  hamiltonian.setStepSize(new_step_size);
  ASSERT_EQ(hamiltonian.getStepSize(), new_step_size);
}

TEST(HamiltonianTests, UpdateDiracUpdatesTheDirac) {
  const Engine engine;

  Hamiltonian hamiltonian(Integrator::LEAPFROG, engine.rng, 0.2);
  auto dirac = DiracOperator(1, 1, 5);
  auto old_dirac_matrix = dirac.getDiracMatrix();
  auto action = Action(-2.7);
  hamiltonian.updateDirac(dirac, action);

  auto new_dirac_matrix = dirac.getDiracMatrix();

  const auto diracs_are_equal = arma::approx_equal(new_dirac_matrix, old_dirac_matrix, "absdiff", 1e-6);

  ASSERT_FALSE(diracs_are_equal) << "The dirac matrix should have been changed";
}