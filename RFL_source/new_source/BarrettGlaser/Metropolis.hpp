//
// Created by Paul Druce on 13/11/2022.
//

#ifndef RFL_METROPOLIS_HPP
#define RFL_METROPOLIS_HPP

#include "Action.hpp"
#include "IAction.hpp"
#include "IAlgorithm.hpp"
#include "IRng.hpp"
#include <armadillo>
#include <memory>

/**
 * Metropolis is a class that encapsulates the implementation of the Metropolis-Hasting
 * for random NCGs.
 *
 */
class Metropolis : public IAlgorithm {
public:
  Metropolis() = delete;
  Metropolis(std::unique_ptr<Action>&& action, double scale, int num_steps, std::unique_ptr<IRng>&& rng)
      : m_action(std::move(action)), m_scale(scale), m_num_steps(num_steps), m_rng(std::move(rng)){};

  void setParams(const double scale, const int iter, std::unique_ptr<IRng>&& rng) {
    this->m_scale = scale;
    this->m_num_steps = iter;
    this->m_rng = std::move(rng);
  }

  double updateDirac(const DiracOperator& dirac) const override {
    return this->run(dirac);
  }

  // MMC routine that performs dual averaging
  double runDualAverage(const DiracOperator& dirac,
                        double target);

protected:
  // MMC routine that doesn't perform dual averaging
  double run(const DiracOperator& dirac) const;

private:
  std::unique_ptr<Action> m_action;
  double m_scale;
  int m_num_steps;
  std::unique_ptr<IRng> m_rng;

  // TODO: This method requires the use of the Action class, not IAction. Refactor this and fix.
  double delta24(const DiracOperator& dirac,
                 const int& x,
                 const int& row_index,
                 const int& column_index,
                 const arma::cx_double& z) const;

  double delta2(const DiracOperator& dirac,
                const int& x,
                const int& row_index,
                const int& column_index,
                const arma::cx_double& z) const;

  double delta4(const DiracOperator& dirac,
                const int& x,
                const int& row_index,
                const int& column_index,
                const arma::cx_double& z) const;

  double runDualAverageCore(const DiracOperator& dirac,
                            const double* s_i,
                            double* s_f) const;

  double runCore(const DiracOperator& dirac,
                 const double* s_i,
                 double* s_f) const;
};

#endif//RFL_METROPOLIS_HPP