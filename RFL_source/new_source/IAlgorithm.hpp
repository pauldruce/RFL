//
// Created by Paul Druce on 31/01/2023.
//

#ifndef RFL_IALGORITHM_HPP
#define RFL_IALGORITHM_HPP

#include "DiracOperator.hpp"
#include "IAction.hpp"

/**
 * This is an Abstract Class/Interface for the various Monte Carlo
 * algorithms that are used in a RFL simulation.
 */
class IAlgorithm {
public:
  /**
   * updateDirac is the only method that an algorithm needs to implement to be used with the RFL
   * library.
   *
   * @param dirac is the DiracOperator the algorithm will operate on.
   * @param action is the Action that the algorithm uses in the MCMC process.
   * @return the acceptance rate of the process.
   */
  virtual double updateDirac(const DiracOperator& dirac, const IAction& action) const = 0;

  virtual ~IAlgorithm() = default;
};

#endif//RFL_IALGORITHM_HPP
