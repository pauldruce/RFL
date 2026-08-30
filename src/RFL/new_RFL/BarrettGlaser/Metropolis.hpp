//
// Created by Paul Druce on 13/11/2022.
//

#ifndef RFL_METROPOLIS_HPP
#define RFL_METROPOLIS_HPP

#include "../IAlgorithm.hpp"
#include "../IRng.hpp"
#include "./Action.hpp"
#include <armadillo>
#include <memory>

/**
 * @class Metropolis
 *
 * @brief Implements the Metropolis algorithm for the Barrett-Glaser action.
 *
 * Samples matrix elements of the Dirac operator according to the probability
 * distribution defined by the Barrett-Glaser action.
 */
class Metropolis final : public IAlgorithm {
public:
  /**
   * Default constructor deleted.
   */
  Metropolis() = delete;

  /**
   * Constructs a Metropolis algorithm instance.
   *
   * @param action Unique pointer to a Barrett-Glaser Action object.
   * @param scale Step scale for matrix element proposals.
   * @param num_steps Number of sweeps per call to updateDirac().
   * @param rng Unique pointer to random number generator engine.
   */
  Metropolis(std::unique_ptr<Action>&& action, const double scale, const int num_steps, std::unique_ptr<IRng>&& rng)
      : m_action(std::move(action)), m_scale(scale), m_num_steps(num_steps), m_rng(std::move(rng)){};

  /**
   * Updates sampling parameters and random number generator.
   *
   * @param scale Step scale for matrix element proposals.
   * @param number_of_steps Number of sweeps per call.
   * @param rng Unique pointer to random number generator engine.
   */
  void setParams(const double scale, const int number_of_steps, std::unique_ptr<IRng>&& rng) {
    this->m_scale = scale;
    this->m_num_steps = number_of_steps;
    this->m_rng = std::move(rng);
  }

  /**
   * Runs Metropolis sweeps on the Dirac operator.
   *
   * @param dirac Dirac operator to update.
   * @return Mean acceptance rate across sweeps.
   */
  double updateDirac(const IDiracOperator& dirac) const override {
    return this->run(dirac);
  }

private:
  std::unique_ptr<Action> m_action;
  double m_scale;
  int m_num_steps;
  std::unique_ptr<IRng> m_rng;

  // MCMC routine without dual-averaging.
  double run(const IDiracOperator& dirac) const;

  // MCMC routine with dual-averaging.
  double runDualAverage(const IDiracOperator& dirac,
                        double target);

  double delta24(const IDiracOperator& dirac,
                 const int& x,
                 const int& row_index,
                 const int& column_index,
                 const arma::cx_double& z) const;

  // TODO: Move to DiracOperator or make static.
  static double delta2(const IDiracOperator& dirac,
                       const int& x,
                       const int& row_index,
                       const int& column_index,
                       const arma::cx_double& z);

  // TODO: Move to DiracOperator or make static.
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
