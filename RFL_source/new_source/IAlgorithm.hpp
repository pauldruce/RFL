//
// Created by Paul Druce on 31/01/2023.
//
/**
 * This is an Abstract Class/Interface for the various Monte Carlo
 * algorithms that are used in a simulation.
 *
 */

#ifndef RFL_IALGORITHM_HPP
#define RFL_IALGORITHM_HPP

#include "Action.hpp"
#include "DiracOperator.hpp"

class IAlgorithm {
public:
  virtual double updateDirac(const DiracOperator& dirac, const Action& action) const = 0;
  virtual ~IAlgorithm() = default;
};

#endif//RFL_IALGORITHM_HPP
