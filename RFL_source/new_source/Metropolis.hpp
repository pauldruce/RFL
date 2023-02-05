//
// Created by Paul Druce on 13/11/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_METROPOLIS_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_METROPOLIS_HPP_

#include "IAlgorithm.hpp"
#include <armadillo>

class Metropolis : public IAlgorithm {
public:
  Metropolis() = default;

  void setParams(const double scale, const int iter, gsl_rng* engine) {
    this->scale = scale;
    this->iter = iter;
    this->engine = engine;
  }

  double updateDirac(const DiracOperator& D, const Action& A) const override {
    return this->MMC(D, A, this->scale, this->iter, this->engine);
  }

private:
  double scale;
  int iter;
  gsl_rng* engine;

  // MMC routine that doesn't perform dual averaging
  double MMC(const DiracOperator& D,
             const Action& A,
             const double& scale,
             const int& iter,
             gsl_rng* engine) const;

  // MMC routine that performs dual averaging
  void MMC_duav(const DiracOperator& D,
                const Action& A,
                double& scale,
                const int& iter,
                gsl_rng* engine,
                const double& target) const;

  double delta24(const DiracOperator& D,
                 const Action& A,
                 const int& x,
                 const int& I,
                 const int& J,
                 const arma::cx_double& z) const;

  double delta2(const DiracOperator& D,
                const int& x,
                const int& I,
                const int& J,
                const arma::cx_double& z) const;

  double delta4(const DiracOperator& D,
                const int& x,
                const int& I,
                const int& J,
                const arma::cx_double& z) const;

  double MMC_duav_core(const DiracOperator& D,
                       const Action& A,
                       const double& scale,
                       gsl_rng* engine,
                       double* s_i,
                       double* s_f) const;

  double MMC_core(const DiracOperator& D,
                  const Action& A,
                  const double& scale,
                  gsl_rng* engine,
                  double* s_i,
                  double* s_f) const;
};

#endif//RFL_RFL_SOURCE_NEW_SOURCE_METROPOLIS_HPP_
