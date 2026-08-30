//
// Created by Paul Druce on 07/12/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_

#include "DiracOperator.hpp"
#include "IAlgorithm.hpp"

/**
 * @class Simulation
 *
 * @brief Manages execution and resources for a Monte Carlo simulation.
 *
 * Combines a Dirac operator and a sampling algorithm to execute simulations.
 */
class Simulation {
public:
  /**
   * Default constructor deleted.
   */
  Simulation() = delete;

  /**
   * Constructs a Simulation instance.
   *
   * @param dirac Unique pointer to a DiracOperator object.
   * @param monte_carlo_algorithm Unique pointer to an IAlgorithm sampling engine.
   */
  Simulation(std::unique_ptr<DiracOperator>&& dirac, std::unique_ptr<IAlgorithm>&& monte_carlo_algorithm);

  /**
   * Runs the simulation algorithm on the Dirac operator.
   *
   * @return Mean acceptance rate across algorithm sweeps or trajectories.
   */
  double run() const {
    return this->m_algorithm->updateDirac(*m_dirac);
  }

  /**
   * Returns a constant reference to the managed Dirac operator.
   *
   * @return Constant reference to DiracOperator.
   */
  const DiracOperator& getDiracOperator() const {
    return *m_dirac;
  }

private:
  std::unique_ptr<DiracOperator> m_dirac;
  std::unique_ptr<IAlgorithm> m_algorithm;
};

#endif//RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
