//
// Created by Paul Druce on 12/11/2022.
//

#include "Hamiltonian.hpp"
#include "../IDiracOperatorDerivatives.hpp"
#include "Action.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

Hamiltonian::Hamiltonian(std::unique_ptr<Action>&& action, const Integrator integrator, const double step_size, std::unique_ptr<IRng>&& rng)
    : m_action(std::move(action)), m_integrator(integrator), m_dt(step_size), m_rng(std::move(rng)) {
}
double Hamiltonian::updateDirac(const IDiracOperator& dirac) const {
  const double acceptance_val_per_iter = this->run(
      dirac,
      10,
      1000);
  return acceptance_val_per_iter;
}

void Hamiltonian::sampleMoments(const IDiracOperator& dirac) const {
  auto& mom = dirac.getMomenta();
  const auto num_matrices = dirac.getNumMatrices();
  const auto mat_dim = dirac.getMatrixDimension();

  for (int i = 0; i < num_matrices; ++i) {
    // loop on indices
    for (int j = 0; j < mat_dim; ++j) {
      const double x = m_rng->getGaussian(1.0);
      mom[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < mat_dim; ++k) {
        const double a = m_rng->getGaussian(1.0);
        const double b = m_rng->getGaussian(1.0);
        mom[i](j, k) = cx_double(a, b) / sqrt(2.);
        mom[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }
}

double Hamiltonian::calculateK(const IDiracOperator& dirac) {
  double res = 0;
  const auto& mom = dirac.getMomenta();
  const auto num_matrices = dirac.getNumMatrices();

  for (int i = 0; i < num_matrices; ++i)
    res += trace(mom[i] * mom[i]).real();

  return res / 2;
}

double Hamiltonian::calculateH(const IDiracOperator& dirac) const {
  return m_action->calculateS(dirac) + calculateK(dirac);
}

void Hamiltonian::leapfrog(const IDiracOperator& dirac,
                           const int& nt,
                           const double g_2) const {
  auto& mat = dirac.getMatrices();
  auto& mom = dirac.getMomenta();
  const auto num_matrices = dirac.getNumMatrices();

  for (int i = 0; i < num_matrices; ++i) {
    mat[i] += m_dt / 2. * mom[i];

    for (int j = 0; j < nt - 1; ++j) {
      mom[i] += -m_dt * derDirac24(dirac, i, true, g_2);
      mat[i] += m_dt * mom[i];
    }

    mom[i] += -m_dt * derDirac24(dirac, i, true, g_2);
    mat[i] += m_dt / 2. * mom[i];
  }
}

void Hamiltonian::omelyan(const IDiracOperator& dirac,
                          const int& nt,
                          const double g_2) const {

  auto& mat = dirac.getMatrices();
  auto& mom = dirac.getMomenta();
  const auto num_matrices = dirac.getNumMatrices();

  for (int i = 0; i < num_matrices; ++i) {
    constexpr double xi = 0.1931833;
    mat[i] += xi * m_dt * mom[i];

    for (int j = 0; j < nt - 1; ++j) {
      mom[i] += -(m_dt / 2.) * derDirac24(dirac, i, true, g_2);
      mat[i] += (1 - 2 * xi) * m_dt * mom[i];
      mom[i] += -(m_dt / 2.) * derDirac24(dirac, i, true, g_2);
      mat[i] += 2 * xi * m_dt * mom[i];
    }

    mom[i] += -(m_dt / 2.) * derDirac24(dirac, i, true, g_2);
    mat[i] += (1 - 2 * xi) * m_dt * mom[i];
    mom[i] += -(m_dt / 2.) * derDirac24(dirac, i, true, g_2);
    mat[i] += xi * m_dt * mom[i];
  }
}

void Hamiltonian::runDualAverage(const IDiracOperator& dirac,
                                 const int& nt,
                                 const int& iter,
                                 const double& target) {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto en_i = vector<double>(4);
  auto en_f = vector<double>(4);

  // dual averaging variables for dt
  double stat = 0;
  const double mu = log(10 * m_dt);
  double log_dt_avg = log(m_dt);

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
    constexpr int i_0 = 10;
    constexpr double kappa = 0.75;
    constexpr double shr = 0.05;
    // if it's not the first iteration set potential to
    // previous final value, otherwise compute it
    if (i) {
      en_i[0] = en_f[0];
      en_i[1] = en_f[1];
    } else {
      en_i[0] = dirac.traceOfDiracSquared();
      en_i[1] = dirac.traceOfDirac4();
    }

    // core part of HMC
    stat += target - runDualAveragingCore(dirac, nt, en_i, en_f);

    // perform dual averaging on dt
    const double log_dt = mu - stat * sqrt(i + 1) / (shr * (i + 1 + i_0));
    m_dt = exp(log_dt);
    const double eta = pow(i + 1, -kappa);
    log_dt_avg = eta * log_dt + (1 - eta) * log_dt_avg;
  }

  // set dt on its final dual averaged value
  m_dt = exp(log_dt_avg);
}

double Hamiltonian::run(const IDiracOperator& dirac,
                        const int& num_iterations,
                        const int& iter) const {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto en_i = vector<double>(4);
  auto en_f = vector<double>(4);

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
      en_i[0] = dirac.traceOfDiracSquared();
      en_i[1] = dirac.traceOfDirac4();
    }

    // core part of HMC
    stat += runCore(dirac, num_iterations, en_i, en_f);
  }

  return stat / iter;
}

double Hamiltonian::run(const IDiracOperator& dirac,
                        const int& nt,
                        const double& dt_min,
                        const double& dt_max,
                        const int& iter) {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto en_i = vector<double>(4);
  auto en_f = vector<double>(4);

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
      en_i[0] = dirac.traceOfDiracSquared();
      en_i[1] = dirac.traceOfDirac4();
    }

    // core part of HMC
    stat += runCore(dirac, nt, dt_min, dt_max, en_i, en_f);
  }

  return stat / iter;
}

double Hamiltonian::runDualAveragingCore(const IDiracOperator& dirac,
                                         const int& nt,
                                         vector<double>& en_i,
                                         vector<double>& en_f) const {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sampleMoments(dirac);

  // store previous configuration
  const auto num_matrices = dirac.getNumMatrices();
  auto mat_bk = vector<cx_mat>(num_matrices);
  auto& mat = dirac.getMatrices();
  for (int j = 0; j < num_matrices; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculateK(dirac);
  en_i[3] = m_action->getG2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (m_integrator == LEAPFROG) {
    leapfrog(dirac, nt, m_action->getG2());
  } else if (m_integrator == OMELYAN) {
    omelyan(dirac, nt, m_action->getG2());
  }

  // calculate final hamiltonian
  en_f[0] = dirac.traceOfDiracSquared();
  en_f[1] = dirac.traceOfDirac4();
  en_f[2] = calculateK(dirac);
  en_f[3] = m_action->getG2() * en_f[0] + en_f[1] + en_f[2];

  // metropolis test

  // sometimes leapfrog diverges and Hf becomes nan.
  // so first of all address this case
  if (std::isnan(en_f[3])) {
    e = 0;
    // restore old configuration
    for (int j = 0; j < num_matrices; ++j)
      mat[j] = mat_bk[j];
    en_f[0] = en_i[0];
    en_f[1] = en_i[1];
    en_f[2] = en_i[2];
    en_f[3] = en_i[3];
  }
  // now do the standard metropolis test
  else if (en_f[3] > en_i[3]) {
    //    double r = gsl_rng_uniform(m_engine);
    const double r = m_rng->getUniform();
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < num_matrices; ++j)
        mat[j] = mat_bk[j];
      en_f[0] = en_i[0];
      en_f[1] = en_i[1];
      en_f[2] = en_i[2];
      en_f[3] = en_i[3];
    }
  }

  return e;
}

double Hamiltonian::runCore(const IDiracOperator& dirac,
                            const int& nt,
                            vector<double>& en_i,
                            vector<double>& en_f) const {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sampleMoments(dirac);

  // store previous configuration
  const auto num_matrices = dirac.getNumMatrices();
  auto mat_bk = vector<cx_mat>(num_matrices);
  auto& mat = dirac.getMatrices();
  for (int j = 0; j < num_matrices; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculateK(dirac);
  en_i[3] = m_action->getG2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (m_integrator == LEAPFROG) {
    leapfrog(dirac, nt, m_action->getG2());
  } else if (m_integrator == OMELYAN) {
    omelyan(dirac, nt, m_action->getG2());
  }

  // calculate final hamiltonian
  en_f[0] = dirac.traceOfDiracSquared();
  en_f[1] = dirac.traceOfDirac4();
  en_f[2] = calculateK(dirac);
  en_f[3] = m_action->getG2() * en_f[0] + en_f[1] + en_f[2];

  // metropolis test
  if (en_f[3] > en_i[3]) {
    const double r = m_rng->getUniform();
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < num_matrices; ++j)
        mat[j] = mat_bk[j];
      en_f[0] = en_i[0];
      en_f[1] = en_i[1];
      en_f[2] = en_i[2];
      en_f[3] = en_i[3];
    }
  }

  return e;
}

double Hamiltonian::runCoreDebug(const IDiracOperator& dirac,
                                 const int& nt) const {
  // exp(-dH) (return value)

  // resample momentum
  sampleMoments(dirac);

  // store previous configuration
  const auto num_matrices = dirac.getNumMatrices();
  auto mat_bk = vector<cx_mat>(num_matrices);
  auto& mat = dirac.getMatrices();
  for (int j = 0; j < num_matrices; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  const double initial_S = m_action->calculateS(dirac);
  const double initial_K = calculateK(dirac);
  const double initial_hamiltonian = initial_S + initial_K;

  // integration
  if (m_integrator == LEAPFROG) {
    leapfrog(dirac, nt, m_action->getG2());
  } else if (m_integrator == OMELYAN) {
    omelyan(dirac, nt, m_action->getG2());
  }

  // calculate final hamiltonian
  const double final_S = m_action->calculateS(dirac);
  const double final_K = calculateK(dirac);
  const double final_hamiltonian = final_S + final_K;

  const double e = exp(initial_hamiltonian - final_hamiltonian);

  // metropolis test
  if (final_hamiltonian > initial_hamiltonian) {
    if (const double r = m_rng->getUniform(); r > e) {
      // restore old configuration
      for (int j = 0; j < num_matrices; ++j)
        mat[j] = mat_bk[j];
    }
  }

  return e;
}

double Hamiltonian::runCore(const IDiracOperator& dirac,
                            const int& nt,
                            const double& dt_min,
                            const double& dt_max,
                            vector<double>& en_i,
                            vector<double>& en_f) {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sampleMoments(dirac);

  // choose uniformly from [dt_min, dt_max)
  this->m_dt = dt_min + (dt_max - dt_min) * m_rng->getUniform();

  // store previous configuration
  const auto num_matrices = dirac.getNumMatrices();
  auto mat_bk = vector<cx_mat>(num_matrices);
  auto& mat = dirac.getMatrices();
  for (int j = 0; j < num_matrices; j++) {
    mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculateK(dirac);
  en_i[3] = m_action->getG2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (m_integrator == LEAPFROG) {
    leapfrog(dirac, nt, m_action->getG2());
  } else if (m_integrator == OMELYAN) {
    omelyan(dirac, nt, m_action->getG2());
  }

  // calculate final hamiltonian
  en_f[0] = dirac.traceOfDiracSquared();
  en_f[1] = dirac.traceOfDirac4();
  en_f[2] = calculateK(dirac);
  en_f[3] = m_action->getG2() * en_f[0] + en_f[1] + en_f[2];

  // metropolis test
  if (en_f[3] > en_i[3]) {
    const double r = m_rng->getUniform();
    e = exp(en_i[3] - en_f[3]);

    if (r > e) {
      // restore old configuration
      for (int j = 0; j < num_matrices; ++j)
        mat[j] = mat_bk[j];
      en_f[0] = en_i[0];
      en_f[1] = en_i[1];
      en_f[2] = en_i[2];
      en_f[3] = en_i[3];
    }
  }

  return e;
}

void Hamiltonian::setStepSize(const double dt) {
  this->m_dt = dt;
}
void Hamiltonian::setIntegrator(const Integrator integrator) { this->m_integrator = integrator; }
