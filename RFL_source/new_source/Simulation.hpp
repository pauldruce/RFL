//
// Created by Paul Druce on 07/12/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_

#include "DiracOperator.hpp"
#include "IAlgorithm.hpp"

/**
 * A class to control the simulation. This class is responsible for the lifespan of
 * of the simulation, when this class goes out of scope, all resources will be cleaned up.
 */
class Simulation {
public:
  /**
   * The default constructor for this class has been disabled.
   */
  Simulation() = delete;

  /**
   * The constructor for the Simulation class. It requires a configured Dirac operator class
   * and a class derived from the IAlgorithm abstract class that implements a Monte Carlo algorithm.
   *
   * @param dirac
   * @param monte_carlo_algorithm
   */
  Simulation(std::unique_ptr<DiracOperator>&& dirac, std::unique_ptr<IAlgorithm>&& monte_carlo_algorithm);

  /**
   * This method starts the simulation and will return when it is complete.
   */
  double run() const {
    return this->m_algorithm->updateDirac(*m_dirac);
  }

  /**
   * This method returns a reference to DiracOperator in use in this class.
   * @return const reference to DiracOperator object
   */
  const DiracOperator& getDiracOperator() const {
    return *m_dirac;
  }

private:
  std::unique_ptr<DiracOperator> m_dirac;
  std::unique_ptr<IAlgorithm> m_algorithm;
};

#endif//RFL_RFL_SOURCE_NEW_SOURCE_SIMULATION_HPP_
