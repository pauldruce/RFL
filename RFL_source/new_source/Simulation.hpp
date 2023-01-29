//
// Created by Paul Druce on 07/12/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#include "DiracOperator.hpp"
#include "Action.hpp"
#include "Metropolis.hpp"

class Simulation {
 public:
  Simulation() = delete;
  Simulation(const DiracOperator &D, const Action &A, const Metropolis &M);
  void run(const double scale, const int stepSize) const;

 private:
  const DiracOperator &m_D;
  const Action &m_A;
  const Metropolis &m_M;
  gsl_rng *engine;

};

#endif //RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
