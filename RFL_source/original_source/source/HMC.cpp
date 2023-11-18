#include "Geom24.hpp"
#include <gsl/gsl_rng.h>

using namespace std;
using namespace arma;

// PD: Moved into HMC class

// HMC routine that performs dual averaging
void Geom24::HMC_duav(const int& Nt,
                      double& dt,
                      const int& iter,
                      gsl_rng* engine,
                      const double& target,
                      const string& integrator) {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto* en_i = new double[4];
  auto* en_f = new double[4];

  // dual averaging variables for dt
  double Stat = 0;
  const double mu = log(10 * dt);
  double log_dt_avg = log(dt);

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
    constexpr int i0 = 10;
    constexpr double kappa = 0.75;
    constexpr double shr = 0.05;
    // if it's not the first interation set potential to
    // previous final value, otherwise compute it
    if (i) {
      en_i[0] = en_f[0];
      en_i[1] = en_f[1];
    } else {
      en_i[0] = dirac2();
      en_i[1] = dirac4();
    }

    // core part of HMC
    Stat += target - HMC_duav_core(Nt, dt, engine, en_i, en_f, integrator);

    // perform dual averaging on dt
    const double log_dt = mu - Stat * sqrt(i + 1) / (shr * (i + 1 + i0));
    dt = exp(log_dt);
    const double eta = pow(i + 1, -kappa);
    log_dt_avg = eta * log_dt + (1 - eta) * log_dt_avg;
  }

  // set dt on its final dual averaged value
  dt = exp(log_dt_avg);

  delete[] en_i;
  delete[] en_f;
}

// HMC routine that doesn't perform dual averaging
double Geom24::HMC(const int& Nt, const double& dt, const int& iter, gsl_rng* engine, const string& integrator) {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto* en_i = new double[4];
  auto* en_f = new double[4];

  // return statistic
  double Stat = 0;

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
    // if it's not the first interation set potential to
    // previous final value, otherwise compute it
    if (i) {
      en_i[0] = en_f[0];
      en_i[1] = en_f[1];
    } else {
      en_i[0] = dirac2();
      en_i[1] = dirac4();
    }

    // core part of HMC
    Stat += HMC_core(Nt, dt, engine, en_i, en_f, integrator);
  }

  delete[] en_i;
  delete[] en_f;

  return (Stat / iter);
}

// HMC routine with randomized integration step
double Geom24::HMC(const int& Nt,
                   const double& dt_min,
                   const double& dt_max,
                   const int& iter,
                   gsl_rng* engine,
                   const string& integrator) {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto* en_i = new double[4];
  auto* en_f = new double[4];

  // return statistic
  double Stat = 0;

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
    // if it's not the first interation set potential to
    // previous final value, otherwise compute it
    if (i) {
      en_i[0] = en_f[0];
      en_i[1] = en_f[1];
    } else {
      en_i[0] = dirac2();
      en_i[1] = dirac4();
    }

    // core part of HMC
    Stat += HMC_core(Nt, dt_min, dt_max, engine, en_i, en_f, integrator);
  }

  delete[] en_i;
  delete[] en_f;

  return (Stat / iter);
}
