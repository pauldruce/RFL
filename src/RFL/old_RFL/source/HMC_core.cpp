#include "Geom24.hpp"
#include <gsl/gsl_rng.h>

using namespace std;
using namespace arma;

double Geom24::HMC_duav_core(const int& Nt,
                             const double& dt,
                             gsl_rng* engine,
                             double* en_i,
                             double* en_f,
                             const string& integrator) {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sample_mom(engine);

  // store previous configuration
  auto* mat_bk = new cx_mat[nHL];
  for (int j = 0; j < nHL; j++)
    mat_bk[j] = mat[j];

  // calculate initial hamiltonian
  en_i[2] = calculate_K();
  en_i[3] = g2 * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (integrator == "leapfrog") {
    leapfrog(Nt, dt);
  } else if (integrator == "omelyan") {
    omelyan(Nt, dt);
  }

  // calculate final hamiltonian
  en_f[0] = dirac2();
  en_f[1] = dirac4();
  en_f[2] = calculate_K();
  en_f[3] = g2 * en_f[0] + en_f[1] + en_f[2];

  // metropolis test

  // sometimes leapfrog diverges and Hf becomes nan.
  // so first of all address this case
  if (std::isnan(en_f[3])) {
    e = 0;
    // restore old configuration
    for (int j = 0; j < nHL; ++j)
      mat[j] = mat_bk[j];
    en_f[0] = en_i[0];
    en_f[1] = en_i[1];
    en_f[2] = en_i[2];
    en_f[3] = en_i[3];
  }
  // now do the standard metropolis test
  else if (en_f[3] > en_i[3]) {
    double r = gsl_rng_uniform(engine);
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < nHL; ++j)
        mat[j] = mat_bk[j];
      en_f[0] = en_i[0];
      en_f[1] = en_i[1];
      en_f[2] = en_i[2];
      en_f[3] = en_i[3];
    }
  }

  delete[] mat_bk;

  return e;
}

double Geom24::HMC_core(const int& Nt,
                        const double& dt,
                        gsl_rng* engine,
                        double* en_i,
                        double* en_f,
                        const string& integrator) {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sample_mom(engine);

  // store previous configuration
  auto* mat_bk = new cx_mat[nHL];
  for (int j = 0; j < nHL; j++)
    mat_bk[j] = mat[j];

  // calculate initial hamiltonian
  en_i[2] = calculate_K();
  en_i[3] = g2 * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (integrator == "leapfrog") {
    leapfrog(Nt, dt);
  } else if (integrator == "omelyan") {
    omelyan(Nt, dt);
  }

  // calculate final hamiltonian
  en_f[0] = dirac2();
  en_f[1] = dirac4();
  en_f[2] = calculate_K();
  en_f[3] = g2 * en_f[0] + en_f[1] + en_f[2];

  // metropolis test
  if (en_f[3] > en_i[3]) {
    double r = gsl_rng_uniform(engine);
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < nHL; ++j)
        mat[j] = mat_bk[j];
      en_f[0] = en_i[0];
      en_f[1] = en_i[1];
      en_f[2] = en_i[2];
      en_f[3] = en_i[3];
    }
  }

  delete[] mat_bk;

  return e;
}

double Geom24::HMC_core_debug(const int& Nt, const double& dt, gsl_rng* engine, const string& integrator) {
  // exp(-dH) (return value)
  double e;

  // resample momentum
  sample_mom(engine);

  // store previous configuration
  auto* mat_bk = new cx_mat[nHL];
  for (int j = 0; j < nHL; j++)
    mat_bk[j] = mat[j];

  // calculate initial hamiltonian
  double Si = calculate_S();
  double Ki = calculate_K();
  double Hi = Si + Ki;

  // integration
  if (integrator == "leapfrog") {
    leapfrog(Nt, dt);
  } else if (integrator == "omelyan") {
    omelyan(Nt, dt);
  }

  // calculate final hamiltonian
  double Sf = calculate_S();
  double Kf = calculate_K();
  double Hf = Sf + Kf;

  e = exp(Hi - Hf);

  // metropolis test
  if (Hf > Hi) {
    double r = gsl_rng_uniform(engine);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < nHL; ++j)
        mat[j] = mat_bk[j];
    }
  }

  delete[] mat_bk;

  return e;
}

double Geom24::HMC_core(const int& Nt,
                        const double& dt_min,
                        const double& dt_max,
                        gsl_rng* engine,
                        double* en_i,
                        double* en_f,
                        const string& integrator) {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sample_mom(engine);

  // choose uniformly from [dt_min, dt_max)
  double dt = dt_min + (dt_max - dt_min) * gsl_rng_uniform(engine);

  // store previous configuration
  auto* mat_bk = new cx_mat[nHL];
  for (int j = 0; j < nHL; j++)
    mat_bk[j] = mat[j];

  // calculate initial hamiltonian
  en_i[2] = calculate_K();
  en_i[3] = g2 * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (integrator == "leapfrog") {
    leapfrog(Nt, dt);
  } else if (integrator == "omelyan") {
    omelyan(Nt, dt);
  }

  // calculate final hamiltonian
  en_f[0] = dirac2();
  en_f[1] = dirac4();
  en_f[2] = calculate_K();
  en_f[3] = g2 * en_f[0] + en_f[1] + en_f[2];

  // metropolis test
  if (en_f[3] > en_i[3]) {
    double r = gsl_rng_uniform(engine);
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < nHL; ++j)
        mat[j] = mat_bk[j];
      en_f[0] = en_i[0];
      en_f[1] = en_i[1];
      en_f[2] = en_i[2];
      en_f[3] = en_i[3];
    }
  }

  delete[] mat_bk;

  return e;
}
