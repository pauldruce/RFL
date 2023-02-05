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
  Simulation(const DiracOperator& D, const Action& A, IAlgorithm& M);
  void run() {
    this->m_M.updateDirac(m_D, m_A);
  }

private:
  const DiracOperator& m_D;
  const Action& m_A;
  IAlgorithm& m_M;
};

#endif//RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
