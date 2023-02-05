//
// Created by Paul Druce on 31/01/2023.
//
/**
 * This is an Abstract Class/Interface for the various Monte Carlo
 * algorithms that are used in a simulation.
 *
 */

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_IALGORITHM_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_IALGORITHM_HPP_

#include "Action.hpp"
#include "DiracOperator.hpp"

class IAlgorithm {
public:
  virtual double updateDirac(const DiracOperator& D, const Action& A) const = 0;
};

#endif//RFL_RFL_SOURCE_NEW_SOURCE_IALGORITHM_HPP_
