//
// Created by Paul Druce on 07/12/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_

#include "DiracOperator.hpp"
#include "IAction.hpp"
#include "IAlgorithm.hpp"
#include <memory>

class Simulation {
public:
  Simulation() = delete;
  Simulation(std::unique_ptr<DiracOperator> &&dirac, std::unique_ptr<IAction> &&action, std::unique_ptr<IAlgorithm> &&monte_carlo_algorithm);
  void run() {
    this->m_algorithm->updateDirac(*m_dirac, *m_action);
  }

private:
  std::unique_ptr<DiracOperator> m_dirac;
  std::unique_ptr<IAction> m_action;
  std::unique_ptr<IAlgorithm> m_algorithm;
};

#endif//RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
