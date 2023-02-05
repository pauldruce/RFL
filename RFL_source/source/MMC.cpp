#include "Geom24.hpp"
#include <cmath>
#include <gsl/gsl_rng.h>
#include <vector>

using namespace std;
using namespace arma;
// PD: All functions moved into Metropolis

// MMC routine that performs dual averaging
void Geom24::MMC_duav(double& scale, const int& iter, gsl_rng* engine, const double& target) {
  // initial (_i) and final (_f) action2 and action4
  auto* s_i = new double[2];
  auto* s_f = new double[2];

  // calculate length of a sweep in terms of dofs
  int Nsw = nHL * dim * dim - nL;

  // dual averaging variables
  const double shr = 0.05;
  const double kappa = 0.75;
  const int i0 = 10;
  double Stat = 0;
  double mu = log(10 * scale);
  double log_scale_avg = log(scale);

  // iter sweeps of metropolis
  for (int i = 0; i < iter; ++i) {
    for (int j = 0; j < Nsw; ++j) {
      // set action to previous final value,
      // unless it's the first iteration
      if (j) {
        s_i[0] = s_f[0];
        s_i[1] = s_f[1];
      } else {
        s_i[0] = dirac2();
        s_i[1] = dirac4();
      }

      Stat += target - MMC_duav_core(scale, engine, s_i, s_f);

      // perform dual averaging
      double log_scale = mu - Stat * sqrt(i + 1) / (shr * (i + 1 + i0));
      scale = exp(log_scale);
      double eta = pow(i + 1, -kappa);
      log_scale_avg = eta * log_scale + (1 - eta) * log_scale_avg;
    }
  }

  // set scale on its final dual averaged value
  scale = exp(log_scale_avg);

  delete[] s_i;
  delete[] s_f;
}

// MMC routine that doesn't perform dual averaging
double Geom24::MMC(const double& scale, const int& iter, gsl_rng* engine) {
  // initial (_i) and final (_f) action2 and action4
  auto* s_i = new double[2];
  auto* s_f = new double[2];

  // calculate length of a sweep in terms of dofs
  int Nsw = nHL * dim * dim - nL;

  // return statistic
  double Stat = 0;

  // iter sweeps of metropolis
  for (int i = 0; i < iter; ++i) {
    for (int j = 0; j < Nsw; ++j) {
      // set action to previous final value,
      // unless it's the first iteration
      if (j) {
        s_i[0] = s_f[0];
        s_i[1] = s_f[1];
      } else {
        s_i[0] = dirac2();
        s_i[1] = dirac4();
      }

      Stat += MMC_core(scale, engine, s_i, s_f);
    }
  }

  delete[] s_i;
  delete[] s_f;

  return (Stat / (iter * Nsw));
}
