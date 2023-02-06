//
// Created by Paul Druce on 07/12/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#include "Action.hpp"
#include "DiracOperator.hpp"
#include "IAlgorithm.hpp"

class Simulation {
public:
  Simulation() = delete;
  Simulation(const DiracOperator& dirac, const Action& action, IAlgorithm& monte_carlo_algorithm);
  void run() {
    this->m_algorithm.updateDirac(m_dirac, m_action);
  }

private:
  const DiracOperator& m_dirac;
  const Action& m_action;
  IAlgorithm& m_algorithm;
};

#endif//RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
