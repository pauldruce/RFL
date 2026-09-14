//
// Created by Paul Druce on 31/01/2023.
//

#ifndef RFL_IALGORITHM_HPP
#define RFL_IALGORITHM_HPP

#include "IDiracOperator.hpp"

/**
 * @interface IAlgorithm
 *
 * @brief Interface for Monte Carlo sampling algorithms in RFL.
 *
 * Defines the contract for updating Dirac operator configurations.
 */
class IAlgorithm {
public:
  /**
   * Updates the Dirac operator and returns the acceptance rate.
   *
   * @param dirac Dirac operator to update.
   * @return Acceptance rate in the range [0, 1].
   */
  virtual double updateDirac(const IDiracOperator& dirac) const = 0;

  virtual ~IAlgorithm() = default;
};

#endif//RFL_IALGORITHM_HPP
