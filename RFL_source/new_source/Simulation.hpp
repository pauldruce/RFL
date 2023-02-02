//
// Created by Paul Druce on 07/12/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#include "DiracOperator.hpp"
#include "Action.hpp"
#include "IAlgorithm.hpp"

class Simulation {
 public:
  Simulation() = delete;
  Simulation(const DiracOperator &D, const Action &A, IAlgorithm &M);
  void run(const double scale, const int stepSize){
	this->m_M.updateDirac(m_D,m_A);
  }

 private:
  const DiracOperator &m_D;
  const Action &m_A;
  IAlgorithm &m_M;
  gsl_rng *engine;

};

#endif //RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
