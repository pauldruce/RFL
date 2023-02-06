//
// Created by Paul Druce on 13/11/2022.
//

#ifndef RFL_METROPOLIS_HPP
#define RFL_METROPOLIS_HPP

#include "IAlgorithm.hpp"
#include <armadillo>

class Metropolis : public IAlgorithm {
public:
  Metropolis() = delete;
  Metropolis(double scale, int iter, gsl_rng* engine)
  : m_scale(scale), m_iter(iter), m_engine(engine){};

  void setParams(const double scale, const int iter, gsl_rng* engine) {
    this->m_scale = scale;
    this->m_iter = iter;
    this->m_engine = engine;
  }

  double updateDirac(const DiracOperator& dirac, const Action& action) const override {
    return this->run(dirac, action, this->m_scale, this->m_iter, this->m_engine);
  }

private:
  double m_scale;
  int m_iter;
  gsl_rng* m_engine;

  // MMC routine that doesn't perform dual averaging
  double run(const DiracOperator& dirac,
             const Action& action,
             const double& scale,
             const int& iter,
             gsl_rng* engine) const;

  // MMC routine that performs dual averaging
  void runDualAverage(const DiracOperator& dirac,
                      const Action& action,
                      double& scale,
                      const int& iter,
                      gsl_rng* engine,
                      const double& target) const;

  double delta24(const DiracOperator& dirac,
                 const Action& action,
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
                            const Action& action,
                            const double& scale,
                            gsl_rng* engine,
                            const double* s_i,
                            double* s_f) const;

  double runCore(const DiracOperator& dirac,
                 const Action& action,
                 const double& scale,
                 gsl_rng* engine,
                 const double* s_i,
                 double* s_f) const;
};

#endif//RFL_METROPOLIS_HPP
