//
// Created by Paul Druce on 12/11/2022.
//

#include "Hamiltonian.hpp"
#include "Action.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

void Hamiltonian::sampleMoments(const DiracOperator& dirac) const {
  auto* mom = dirac.get_moms();
  for (int i = 0; i < dirac.nHL; ++i) {
    // loop on indices
    for (int j = 0; j < dirac.dim; ++j) {
      double x;
      x = gsl_ran_gaussian(m_engine, 1.);
      mom[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < dirac.dim; ++k) {
        double a, b;
        a = gsl_ran_gaussian(m_engine, 1.);
        b = gsl_ran_gaussian(m_engine, 1.);
        mom[i](j, k) = cx_double(a, b) / sqrt(2.);
        mom[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }
}

double Hamiltonian::calculateK(const DiracOperator& dirac) const {
  double res = 0;
  auto* mom = dirac.get_moms();

  for (int i = 0; i < dirac.nHL; ++i)
    res += trace(mom[i] * mom[i]).real();

  return res / 2;
}

double Hamiltonian::calculateH(const DiracOperator& dirac, const Action& action) const {
  return action.calculateS(dirac) + calculateK(dirac);
}

void Hamiltonian::leapfrog(const DiracOperator& dirac,
                           const int& nt,
                           double g_2) const {
  auto* mat = dirac.get_mats();
  auto* mom = dirac.get_moms();

  for (int i = 0; i < dirac.nHL; ++i) {
    mat[i] += (m_dt / 2.) * mom[i];

    for (int j = 0; j < nt - 1; ++j) {
      mom[i] += -m_dt * dirac.der_dirac24(i, true, g_2);
      mat[i] += m_dt * mom[i];
    }

    mom[i] += -m_dt * dirac.der_dirac24(i, true, g_2);
    mat[i] += (m_dt / 2.) * mom[i];
  }
}

void Hamiltonian::omelyan(const DiracOperator& dirac,
                          const int& nt,
                          double g_2) const {
  double xi = 0.1931833;

  auto* mat = dirac.get_mats();
  auto* mom = dirac.get_moms();

  for (int i = 0; i < dirac.nHL; ++i) {
    mat[i] += xi * m_dt * mom[i];

    for (int j = 0; j < nt - 1; ++j) {
      mom[i] += -(m_dt / 2.) * dirac.der_dirac24(i, true, g_2);
      mat[i] += (1 - 2 * xi) * m_dt * mom[i];
      mom[i] += -(m_dt / 2.) * dirac.der_dirac24(i, true, g_2);
      mat[i] += 2 * xi * m_dt * mom[i];
    }

    mom[i] += -(m_dt / 2.) * dirac.der_dirac24(i, true, g_2);
    mat[i] += (1 - 2 * xi) * m_dt * mom[i];
    mom[i] += -(m_dt / 2.) * dirac.der_dirac24(i, true, g_2);
    mat[i] += xi * m_dt * mom[i];
  }
}

void Hamiltonian::runHmcDuav(const DiracOperator& dirac,
                             const Action& action,
                             const int& nt,
                             const int& iter,
                             const double& target) {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto* en_i = new double[4];
  auto* en_f = new double[4];

  // dual averaging variables for dt
  const double shr = 0.05;
  const double kappa = 0.75;
  const int i_0 = 10;
  double stat = 0;
  double mu = log(10 * m_dt);
  double log_dt_avg = log(m_dt);

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
    // if it's not the first iteration set potential to
    // previous final value, otherwise compute it
    if (i) {
      en_i[0] = en_f[0];
      en_i[1] = en_f[1];
    } else {
      en_i[0] = action.dirac2(dirac);
      en_i[1] = action.dirac4(dirac);
    }

    // core part of HMC
    stat += target - runHmcDuavCore(dirac, action, nt, en_i, en_f);

    // perform dual averaging on dt
    double log_dt = mu - stat * sqrt(i + 1) / (shr * (i + 1 + i_0));
    m_dt = exp(log_dt);
    double eta = pow(i + 1, -kappa);
    log_dt_avg = eta * log_dt + (1 - eta) * log_dt_avg;
  }

  // set dt on its final dual averaged value
  m_dt = exp(log_dt_avg);

  delete[] en_i;
  delete[] en_f;
}

double Hamiltonian::runHmc(const DiracOperator& dirac,
                           const Action& action,
                           const int& num_iterations,
                           const int& iter) const {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto* en_i = new double[4];
  auto* en_f = new double[4];

  // return statistic
  double stat = 0;

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
    // if it's not the first iteration set potential to
    // previous final value, otherwise compute it
    if (i) {
      en_i[0] = en_f[0];
      en_i[1] = en_f[1];
    } else {
      en_i[0] = action.dirac2(dirac);
      en_i[1] = action.dirac4(dirac);
    }

    // core part of HMC
    stat += runHmcCore(dirac, action, num_iterations, en_i, en_f);
  }

  delete[] en_i;
  delete[] en_f;

  return (stat / iter);
}

double Hamiltonian::runHmc(const DiracOperator& dirac,
                           const Action& action,
                           const int& nt,
                           const double& dt_min,
                           const double& dt_max,
                           const int& iter) {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto* en_i = new double[4];
  auto* en_f = new double[4];

  // return statistic
  double stat = 0;

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
    // if it's not the first iteration set potential to
    // previous final value, otherwise compute it
    if (i) {
      en_i[0] = en_f[0];
      en_i[1] = en_f[1];
    } else {
      en_i[0] = action.dirac2(dirac);
      en_i[1] = action.dirac4(dirac);
    }

    // core part of HMC
    stat += runHmcCore(dirac, action, nt, dt_min, dt_max, en_i, en_f);
  }

  delete[] en_i;
  delete[] en_f;

  return (stat / iter);
}

double Hamiltonian::runHmcDuavCore(const DiracOperator& dirac,
                                   const Action& action,
                                   const int& nt,
                                   double* en_i,
                                   double* en_f) const {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sampleMoments(dirac);

  // store previous configuration
  auto* mat_bk = new cx_mat[dirac.nHL];
  auto* mat = dirac.get_mats();
  for (int j = 0; j < dirac.nHL; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculateK(dirac);
  en_i[3] = action.getG2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (m_integrator == Integrator::LEAPFROG) {
    leapfrog(dirac, nt, action.getG2());
  } else if (m_integrator == Integrator::OMELYAN) {
    omelyan(dirac, nt, action.getG2());
  }

  // calculate final hamiltonian
  en_f[0] = action.dirac2(dirac);
  en_f[1] = action.dirac4(dirac);
  en_f[2] = calculateK(dirac);
  en_f[3] = action.getG2() * en_f[0] + en_f[1] + en_f[2];

  // metropolis test

  // sometimes leapfrog diverges and Hf becomes nan.
  // so first of all address this case
  if (std::isnan(en_f[3])) {
    e = 0;
    // restore old configuration
    for (int j = 0; j < dirac.nHL; ++j)
      mat[j] = mat_bk[j];
    en_f[0] = en_i[0];
    en_f[1] = en_i[1];
    en_f[2] = en_i[2];
    en_f[3] = en_i[3];
  }
    // now do the standard metropolis test
  else if (en_f[3] > en_i[3]) {
    double r = gsl_rng_uniform(m_engine);
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < dirac.nHL; ++j)
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

double Hamiltonian::runHmcCore(const DiracOperator& dirac,
                               const Action& action,
                               const int& nt,
                               double* en_i,
                               double* en_f) const {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sampleMoments(dirac);

  // store previous configuration
  auto* mat_bk = new cx_mat[dirac.nHL];
  auto* mat = dirac.get_mats();
  for (int j = 0; j < dirac.nHL; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculateK(dirac);
  en_i[3] = action.getG2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (m_integrator == Integrator::LEAPFROG) {
    leapfrog(dirac, nt, action.getG2());
  } else if (m_integrator == Integrator::OMELYAN) {
    omelyan(dirac, nt, action.getG2());
  }

  // calculate final hamiltonian
  en_f[0] = action.dirac2(dirac);
  en_f[1] = action.dirac4(dirac);
  en_f[2] = calculateK(dirac);
  en_f[3] = action.getG2() * en_f[0] + en_f[1] + en_f[2];

  // metropolis test
  if (en_f[3] > en_i[3]) {
    double r = gsl_rng_uniform(m_engine);
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < dirac.nHL; ++j)
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

double Hamiltonian::runHmcCoreDebug(const DiracOperator& dirac,
                                    const Action& action,
                                    const int& nt) const {
  // exp(-dH) (return value)
  double e;

  // resample momentum
  sampleMoments(dirac);

  // store previous configuration
  auto* mat_bk = new cx_mat[dirac.nHL];
  auto* mat = dirac.get_mats();
  for (int j = 0; j < dirac.nHL; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  double initial_S = action.calculateS(dirac);
  double initial_K = calculateK(dirac);
  double initial_hamiltonian = initial_S + initial_K;

  // integration
  if (m_integrator == Integrator::LEAPFROG) {
    leapfrog(dirac, nt, action.getG2());
  } else if (m_integrator == Integrator::OMELYAN) {
    omelyan(dirac, nt, action.getG2());
  }

  // calculate final hamiltonian
  double final_S = action.calculateS(dirac);
  double final_K = calculateK(dirac);
  double final_hamiltonian = final_S + final_K;

  e = exp(initial_hamiltonian - final_hamiltonian);

  // metropolis test
  if (final_hamiltonian > initial_hamiltonian) {
    double r = gsl_rng_uniform(m_engine);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < dirac.nHL; ++j)
        mat[j] = mat_bk[j];
    }
  }

  delete[] mat_bk;

  return e;
}

double Hamiltonian::runHmcCore(const DiracOperator& dirac,
                               const Action& action,
                               const int& nt,
                               const double& dt_min,
                               const double& dt_max,
                               double* en_i,
                               double* en_f) {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sampleMoments(dirac);

  // choose uniformly from [dt_min, dt_max)
  this->m_dt = dt_min + (dt_max - dt_min) * gsl_rng_uniform(m_engine);

  // store previous configuration
  auto* mat_bk = new cx_mat[dirac.nHL];
  auto* mat = dirac.get_mats();
  for (int j = 0; j < dirac.nHL; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculateK(dirac);
  en_i[3] = action.getG2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (m_integrator == Integrator::LEAPFROG) {
    leapfrog(dirac, nt, action.getG2());
  } else if (m_integrator == Integrator::OMELYAN) {
    omelyan(dirac, nt, action.getG2());
  }

  // calculate final hamiltonian
  en_f[0] = action.dirac2(dirac);
  en_f[1] = action.dirac4(dirac);
  en_f[2] = calculateK(dirac);
  en_f[3] = action.getG2() * en_f[0] + en_f[1] + en_f[2];

  // metropolis test
  if (en_f[3] > en_i[3]) {
    double r = gsl_rng_uniform(m_engine);
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < dirac.nHL; ++j)
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

void Hamiltonian::setStepSize(double dt) {
  this->m_dt = dt;
}
void Hamiltonian::setIntegrator(Integrator integrator) { this->m_integrator = integrator; }
void Hamiltonian::setEngine(const gsl_rng* engine) {
  this->m_engine = engine;
}

Hamiltonian::Hamiltonian(Integrator integrator, const gsl_rng* engine, double step_size)
    : m_integrator(integrator), m_engine(engine), m_dt(step_size) {
}
double Hamiltonian::updateDirac(const DiracOperator& dirac, const Action& action) const {
  const double acceptance_val_per_iter = this->runHmc(
      dirac, action,
      10,
      1000);
  return acceptance_val_per_iter;
}