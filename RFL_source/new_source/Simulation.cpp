//
// Created by Paul Druce on 07/12/2022.
//

#include "Simulation.hpp"

Simulation::Simulation(const DiracOperator &D, const Action &A, IAlgorithm &M) :
	m_D(D), m_A(A), m_M(M) {
  // Initialize the random number generator
  this->engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));
}
