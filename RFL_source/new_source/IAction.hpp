//
// Created by Paul Druce on 15/03/2023.
//

#ifndef RFL_IACTION_HPP
#define RFL_IACTION_HPP

#include "IDiracOperator.hpp"

/**
 * IAction is an abstract class that provides an interface for any action you may want to use in a
 * random NCG simulation.
 */
class IAction {
public:
  /**
    * calculateS is a pure virtual method that needs to be implemented by a derived class
    *
    * Any implementation should calculate the action and return its value as a double precision
    * floating point number.
    *
    * @param dirac is the Dirac operator that the action should be calculated for.
    * @return the value of the action for the input Dirac operator.
    */
  virtual double calculateS(const IDiracOperator& dirac) const = 0;

  virtual ~IAction() = default;
};

#endif//RFL_IACTION_HPP
