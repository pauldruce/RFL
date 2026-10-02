//
// Created by Paul Druce on 13/11/2022.
//
#include "./Metropolis.hpp"
#include <cmath>

using namespace std;
using namespace arma;

double Metropolis::runDualAverage(
    IDiracOperator& dirac,
    const double target) {
  // Initial and final action traces.
  double s_i[2]{0.0, 0.0};
  double s_f[2]{0.0, 0.0};
  const auto mat_dim = dirac.getMatrixDimension();
  // Number of degrees of freedom per sweep:
  // Each matrix has mat_dim^2 elements; anti-Hermitian matrices have constrained diagonal degrees of freedom.
  const int nsw = dirac.getNumMatrices() * mat_dim * mat_dim - dirac.getNumAntiHermitianMatrices();

  // Dual-averaging variables.
  double stat = 0;
  const double mu = log(10 * m_scale);
  double log_scale_avg = log(m_scale);

  // Run sweeps of the Metropolis algorithm.
  for (int i = 0; i < m_num_steps; ++i) {
    for (int j = 0; j < nsw; ++j) {
      constexpr int i_0 = 10;
      constexpr double kappa = 0.75;
      constexpr double shr = 0.05;
      // Set action to previous final value, unless on the first step of the sweep.
      if (j) {
        s_i[0] = s_f[0];
        s_i[1] = s_f[1];
      } else {
        s_i[0] = dirac.traceOfDiracSquared();
        s_i[1] = dirac.traceOfDirac4();
      }

      stat += target - runDualAverageCore(dirac, s_i, s_f);

      // Perform dual-averaging.
      const double log_scale = mu - stat * sqrt(i + 1) / (shr * (i + 1 + i_0));
      m_scale = exp(log_scale);
      const double eta = pow(i + 1, -kappa);
      log_scale_avg = eta * log_scale + (1 - eta) * log_scale_avg;
    }
  }

  // Set proposal scale to its final dual-averaged value.
  m_scale = exp(log_scale_avg);

  return (stat / (m_num_steps * nsw));
}

double Metropolis::run(IDiracOperator& dirac) const {
  // Initial and final action traces.
  double s_i[2]{0.0, 0.0};
  double s_f[2]{0.0, 0.0};
  const auto mat_dim = dirac.getMatrixDimension();
  // Number of degrees of freedom per sweep:
  // Each matrix has mat_dim^2 elements; anti-Hermitian matrices have constrained diagonal degrees of freedom.
  const int nsw = dirac.getNumMatrices() * mat_dim * mat_dim - dirac.getNumAntiHermitianMatrices();

  // Acceptance statistic accumulator.
  double stat = 0;

  // Run sweeps of the Metropolis algorithm.
  for (int i = 0; i < m_num_steps; ++i) {
    for (int j = 0; j < nsw; ++j) {
      // Set action to previous final value, unless on the first step of the sweep.
      if (j) {
        s_i[0] = s_f[0];
        s_i[1] = s_f[1];
      } else {
        s_i[0] = dirac.traceOfDiracSquared();
        s_i[1] = dirac.traceOfDirac4();
      }

      stat += runCore(dirac, s_i, s_f);
    }
  }

  return (stat / (m_num_steps * nsw));
}

double Metropolis::runDualAverageCore(
    IDiracOperator& dirac,
    const double* s_i,
    double* s_f) const {
  // Acceptance probability.
  double e;
  const auto num_matrices = dirac.getNumMatrices();
  const auto mat_dim = dirac.getMatrixDimension();

  // Propose matrix element update:
  // Off-diagonal perturbations are complex: delta_M = z * E_{ij} + conj(z) * E_{ji}.
  // Diagonal perturbations are purely real (Im(z) = 0) to preserve Hermiticity.
  const int x = static_cast<int>(m_rng->getUniformInt(0, num_matrices - 1));
  const int row_index = static_cast<int>(m_rng->getUniformInt(0, mat_dim - 1));
  const int column_index = static_cast<int>(m_rng->getUniformInt(0, mat_dim - 1));

  double re = 0;
  cx_double z;
  if (row_index != column_index) {
    double im = 0;
    re = m_scale * (-1. + 2. * m_rng->getUniform());
    im = m_scale * (-1. + 2. * m_rng->getUniform());
    z = cx_double(re, im);
  } else {
    re = m_scale * (-1. + 2. * m_rng->getUniform());
    z = cx_double(re, 0);
  }

  const double delta_2 = m_action->delta2(dirac, x, row_index, column_index, z);
  const double delta_4 = m_action->delta4(dirac, x, row_index, column_index, z);
  const double action_delta = m_action->delta24(delta_2, delta_4);

  auto& mat = dirac.getMatrices();
  // Metropolis accept/reject test.
  if (action_delta < 0) {
    // For off-diagonal moves, delta_M = z * E_{ij} + conj(z) * E_{ji}.
    // For diagonal moves, delta_M = z + conj(z) = 2 * Re(z). Using 2. * z.real() guarantees
    // the diagonal remains strictly real.
    if (row_index != column_index) {
      mat[x](row_index, column_index) += z;
      mat[x](column_index, row_index) += conj(z);
    } else {
      mat[x](row_index, row_index) += 2. * z.real();
    }

    // Update action.
    s_f[0] = s_i[0] + delta_2;
    s_f[1] = s_i[1] + delta_4;

    // Proposal accepted.
    e = 1;
  } else {
    e = exp(-action_delta);
    const double p = m_rng->getUniform();

    if (e > p) {
      // For off-diagonal moves, delta_M = z * E_{ij} + conj(z) * E_{ji}.
      // For diagonal moves, delta_M = z + conj(z) = 2 * Re(z). Using 2. * z.real() guarantees
      // the diagonal remains strictly real.
      if (row_index != column_index) {
        mat[x](row_index, column_index) += z;
        mat[x](column_index, row_index) += conj(z);
      } else {
        mat[x](row_index, row_index) += 2. * z.real();
      }

      // Update action.
      s_f[0] = s_i[0] + delta_2;
      s_f[1] = s_i[1] + delta_4;
    } else {
      s_f[0] = s_i[0];
      s_f[1] = s_i[1];
    }
  }

  return e;
}

double Metropolis::runCore(
    IDiracOperator& dirac,
    const double* s_i,
    double* s_f) const {
  // Acceptance probability.
  double ret = 0;

  const auto num_matrices = dirac.getNumMatrices();
  const auto mat_dim = dirac.getMatrixDimension();

  // Propose matrix element update:
  // Off-diagonal perturbations are complex: delta_M = z * E_{ij} + conj(z) * E_{ji}.
  // Diagonal perturbations are purely real (Im(z) = 0) to preserve Hermiticity.
  const int x = static_cast<int>(m_rng->getUniformInt(0, num_matrices - 1));
  const int row_index = static_cast<int>(m_rng->getUniformInt(0, mat_dim - 1));
  const int column_index = static_cast<int>(m_rng->getUniformInt(0, mat_dim - 1));
  double re = 0;
  cx_double z;
  if (row_index != column_index) {
    double im = 0;
    re = m_scale * (-1. + 2. * m_rng->getUniform());
    im = m_scale * (-1. + 2. * m_rng->getUniform());
    z = cx_double(re, im);
  } else {
    re = m_scale * (-1. + 2. * m_rng->getUniform());
    z = cx_double(re, 0);
  }

  const double delta_2 = m_action->delta2(dirac, x, row_index, column_index, z);
  const double delta_4 = m_action->delta4(dirac, x, row_index, column_index, z);
  const double action_delta = m_action->delta24(delta_2, delta_4);

  auto& mat = dirac.getMatrices();
  // Metropolis accept/reject test.
  if (action_delta < 0) {
    // For off-diagonal moves, delta_M = z * E_{ij} + conj(z) * E_{ji}.
    // For diagonal moves, delta_M = z + conj(z) = 2 * Re(z). Using 2. * z.real() guarantees
    // the diagonal remains strictly real.
    if (row_index != column_index) {
      mat[x](row_index, column_index) += z;
      mat[x](column_index, row_index) += conj(z);
    } else {
      mat[x](row_index, row_index) += 2. * z.real();
    }

    // Update action.
    s_f[0] = s_i[0] + delta_2;
    s_f[1] = s_i[1] + delta_4;

    // Proposal accepted.
    ret = 1;
  } else {
    const double e = exp(-action_delta);
    const double p = m_rng->getUniform();

    if (e > p) {
      // For off-diagonal moves, delta_M = z * E_{ij} + conj(z) * E_{ji}.
      // For diagonal moves, delta_M = z + conj(z) = 2 * Re(z). Using 2. * z.real() guarantees
      // the diagonal remains strictly real.
      if (row_index != column_index) {
        mat[x](row_index, column_index) += z;
        mat[x](column_index, row_index) += conj(z);
      } else {
        mat[x](row_index, row_index) += 2. * z.real();
      }

      // Update action.
      s_f[0] = s_i[0] + delta_2;
      s_f[1] = s_i[1] + delta_4;

      // Proposal accepted.
      ret = 1;
    } else {
      s_f[0] = s_i[0];
      s_f[1] = s_i[1];
    }
  }

  return ret;
}
