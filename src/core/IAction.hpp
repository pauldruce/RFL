//
// Created by Paul Druce on 15/03/2023.
//

#ifndef RFL_IACTION_HPP
#define RFL_IACTION_HPP

#include "IDiracOperator.hpp"

/**
 * @interface IAction
 *
 * @brief Interface for spectral action functionals in RFL.
 *
 * Defines the contract for evaluating an action functional on a Dirac operator.
 */
class IAction {
public:
  /**
   * Calculates the action value for the given Dirac operator.
   *
   * @param dirac Dirac operator reference.
   * @return Action value as a double.
   */
  virtual double calculateS(const IDiracOperator& dirac) const = 0;

  virtual ~IAction() = default;
};

#endif//RFL_IACTION_HPP
