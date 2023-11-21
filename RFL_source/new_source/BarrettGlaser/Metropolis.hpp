//
// Created by Paul Druce on 13/11/2022.
//

#ifndef RFL_METROPOLIS_HPP
#define RFL_METROPOLIS_HPP

#include "Action.hpp"
#include "IAlgorithm.hpp"
#include "IRng.hpp"
#include <armadillo>
#include <memory>

/**
 * Metropolis is a class that encapsulates the implementation of the Metropolis-Hasting
 * algorithm for random non-commutative geometries that are following a simulation
 * governed by the Barrett-Glaser action.
 */
class Metropolis final : public IAlgorithm {
public:
  /**
   * The default constructor has been disabled, please use another constructor.
   */
  Metropolis() = delete;

  /**
   * @param action is a unique_ptr to a Barrett-Glaser action class, Action.
   * @param scale is a parameter to control to step size taken in the Metropolis algorithm
   * @param num_steps is the number of steps to take for each call to the updateDirac method below.
   * @param rng is a unique_ptr to an implementation of the abstract rng class IRng.
   */
  Metropolis(std::unique_ptr<Action>&& action, const double scale, const int num_steps, std::unique_ptr<IRng>&& rng)
      : m_action(std::move(action)), m_scale(scale), m_num_steps(num_steps), m_rng(std::move(rng)){};

  /**
   * setParams updates the scale, number_of_steps and rng parameters that are passed in as part of
   * the constructor.
   */
  void setParams(const double scale, const int number_of_steps, std::unique_ptr<IRng>&& rng) {
    this->m_scale = scale;
    this->m_num_steps = number_of_steps;
    this->m_rng = std::move(rng);
  }

  double updateDirac(const IDiracOperator& dirac) const override {
    return this->run(dirac);
  }

private:
  std::unique_ptr<Action> m_action;
  double m_scale;
  int m_num_steps;
  std::unique_ptr<IRng> m_rng;

  // MMC routine that doesn't perform dual averaging
  double run(const IDiracOperator& dirac) const;

  // MMC routine that performs dual averaging
  double runDualAverage(const IDiracOperator& dirac,
                        double target);

  double delta24(const IDiracOperator& dirac,
                 const int& x,
                 const int& row_index,
                 const int& column_index,
                 const arma::cx_double& z) const;

  // TODO: These can be made static or members of DiracOperator
  static double delta2(const IDiracOperator& dirac,
                       const int& x,
                       const int& row_index,
                       const int& column_index,
                       const arma::cx_double& z);

  // TODO: These can be made static or members of DiracOperator
  static double delta4(const IDiracOperator& dirac,
                       const int& x,
                       const int& row_index,
                       const int& column_index,
                       const arma::cx_double& z);

  double runDualAverageCore(const IDiracOperator& dirac,
                            const double* s_i,
                            double* s_f) const;

  double runCore(const IDiracOperator& dirac,
                 const double* s_i,
                 double* s_f) const;
};

#endif//RFL_METROPOLIS_HPP
