//
// Created by Paul Druce on 13/11/2022.
//
#include "Metropolis.hpp"
#include <cmath>

using namespace std;
using namespace arma;

double Metropolis::delta24(const IDiracOperator& dirac,
                           const int& x,
                           const int& row_index,
                           const int& column_index,
                           const cx_double& z) const {
  return m_action->getG2() * delta2(dirac, x, row_index, column_index, z) + delta4(dirac, x, row_index, column_index, z);
}

double Metropolis::delta2(const IDiracOperator& dirac,
                          const int& x,
                          const int& row_index,
                          const int& column_index,
                          const cx_double& z) {
  const auto& mat = dirac.getMatrices();
  const auto& eps = dirac.getEpsilons();
  const auto mat_dim = dirac.getMatrixDimension();
  const auto gamma_dim = dirac.getGammaDimension();

  if (row_index != column_index) {
    return 4. * gamma_dim * mat_dim * (2. * (z * mat[x](column_index, row_index)).real() + norm(z));
  } else {
    const double tr_m = trace(mat[x]).real();
    return 8. * gamma_dim * z.real() * (mat_dim * (mat[x](row_index, row_index).real() + z.real()) + eps[x] * (tr_m + z.real()));
  }
}

double Metropolis::delta4(const IDiracOperator& dirac,
                          const int& x,
                          const int& row_index,
                          const int& column_index,
                          const cx_double& z) {
  double res = 0.;

  auto& omega_table_4 = dirac.getOmegaTable4();
  auto& mat = dirac.getMatrices();
  auto& eps = dirac.getEpsilons();
  auto mat_dim = dirac.getMatrixDimension();
  auto gamma_dim = dirac.getGammaDimension();
  auto num_matrices = dirac.getNumMatrices();

  // D^3 dD part
  for (int i_3 = 0; i_3 < num_matrices; ++i_3) {
    for (int i_2 = 0; i_2 < num_matrices; ++i_2) {
      for (int i_1 = 0; i_1 <= i_3; ++i_1) {
        cx_double cliff = omega_table_4[x + num_matrices * (i_3 + num_matrices * (i_2 + num_matrices * i_1))];

        if (fabs(cliff.real()) > 1e-10 || fabs(cliff.imag()) > 1e-10) {
          // compute necessary matrix products
          cx_mat m_1_m_2 = mat[i_1] * mat[i_2];
          cx_mat m_2_m_3 = mat[i_2] * mat[i_3];
          cx_mat m_1_m_3 = mat[i_1] * mat[i_3];
          cx_mat m_1_m_2_m_3 = mat[i_1] * m_2_m_3;

          // compute necessary traces
          double tr_m_1 = trace(mat[i_1]).real();
          double tr_m_2 = trace(mat[i_2]).real();
          double tr_m_3 = trace(mat[i_3]).real();
          double tr_m_1_m_2 = trace(m_1_m_2).real();
          double tr_m_2_m_3 = trace(m_2_m_3).real();
          double tr_m_1_m_3 = trace(m_1_m_3).real();
          cx_double tr_m_1_m_2_m_3 = trace(m_1_m_2_m_3);

          // off-diagonal update
          if (row_index != column_index) {

            // compute terms
            // _______________________________________________________________________________________
            cx_double t_1 = m_1_m_2_m_3(column_index, row_index) * z + m_1_m_2_m_3(row_index, column_index) * conj(z);
            t_1 = t_1 + conj(t_1) * (double)(eps[i_1] * eps[i_2] * eps[i_3] * eps[x]);
            t_1 *= (double)mat_dim;

            cx_double t_2 = m_1_m_2(column_index, row_index) * z + m_1_m_2(row_index, column_index) * conj(z);
            t_2 = t_2 * (double)(eps[i_3]) + conj(t_2) * (double)(eps[i_1] * eps[i_2] * eps[x]);
            t_2 = t_2 * tr_m_3;
            t_1 += t_2;

            cx_double t_3 = m_1_m_3(column_index, row_index) * z + m_1_m_3(row_index, column_index) * conj(z);
            t_3 = t_3 * (double)(eps[i_2]) + conj(t_3) * (double)(eps[i_1] * eps[i_3] * eps[x]);
            t_3 = t_3 * tr_m_2;
            t_1 += t_3;

            cx_double t_4 = m_2_m_3(column_index, row_index) * z + m_2_m_3(row_index, column_index) * conj(z);
            t_4 = t_4 * (double)(eps[i_1]) + conj(t_4) * (double)(eps[i_2] * eps[i_3] * eps[x]);
            t_4 = t_4 * tr_m_1;
            t_1 += t_4;

            double t_5 = tr_m_1_m_2 * (eps[i_1] * eps[i_2] + eps[i_3] * eps[x]);
            t_5 *= 2. * (mat[i_3](column_index, row_index) * z).real();
            t_1 += t_5;

            double t_6 = tr_m_2_m_3 * (eps[i_2] * eps[i_3] + eps[i_1] * eps[x]);
            t_6 *= 2. * (mat[i_1](column_index, row_index) * z).real();
            t_1 += t_6;

            double t_7 = tr_m_1_m_3 * (eps[i_1] * eps[i_3] + eps[i_2] * eps[x]);
            t_7 *= 2. * (mat[i_2](column_index, row_index) * z).real();
            t_1 += t_7;
            //________________________________________________________________________________________

            // add to total
            if (i_1 != i_3) {
              res += 2. * (cliff * t_1).real();
            } else {
              res += (cliff * t_1).real();
            }
          }

          // diagonal update
          else {

            // compute terms
            // _______________________________________________________________________________________
            cx_double t_1 = m_1_m_2_m_3(row_index, row_index);
            t_1 = t_1 + conj(t_1) * (double)(eps[i_1] * eps[i_2] * eps[i_3] * eps[x]);
            t_1 = t_1 * (double)mat_dim;

            cx_double t_2 = m_1_m_2(row_index, row_index);
            t_2 = t_2 * (double)(eps[i_3]) + conj(t_2) * (double)(eps[i_1] * eps[i_2] * eps[x]);
            t_2 *= tr_m_3;
            t_1 += t_2;

            cx_double t_3 = m_1_m_3(row_index, row_index);
            t_3 = t_3 * (double)(eps[i_2]) + conj(t_3) * (double)(eps[i_1] * eps[i_3] * eps[x]);
            t_3 *= tr_m_2;
            t_1 += t_3;

            cx_double t_4 = m_2_m_3(row_index, row_index);
            t_4 = t_4 * (double)(eps[i_1]) + conj(t_4) * (double)(eps[i_2] * eps[i_3] * eps[x]);
            t_4 *= tr_m_1;
            t_1 += t_4;

            double t_5 = tr_m_1_m_2 * (eps[i_1] * eps[i_2] + eps[i_3] * eps[x]);
            t_5 *= mat[i_3](row_index, row_index).real();
            t_1 += t_5;

            double t_6 = tr_m_2_m_3 * (eps[i_2] * eps[i_3] + eps[i_1] * eps[x]);
            t_6 *= mat[i_1](row_index, row_index).real();
            t_1 += t_6;

            double t_7 = tr_m_1_m_3 * (eps[i_1] * eps[i_3] + eps[i_2] * eps[x]);
            t_7 *= mat[i_2](row_index, row_index).real();
            t_1 += t_7;

            cx_double t_8 =
                conj(tr_m_1_m_2_m_3) * (double)(eps[i_1] * eps[i_2] * eps[i_3]) + tr_m_1_m_2_m_3 * (double)(eps[x]);
            t_1 += t_8;
            //________________________________________________________________________________________

            // add to total
            if (i_1 != i_3) {
              res += (cliff * t_1).real() * 4. * z.real();
            } else {
              res += (cliff * t_1).real() * 2. * z.real();
            }
          }
        }
      }
    }
  }

  res *= 4.;

  // D^2 dD^2 and D dD D dD term
  double temp = 0;
  for (int i = 0; i < num_matrices; ++i) {
    double cliff = omega_table_4[x + num_matrices * (i + num_matrices * (x + num_matrices * i))].real();

    // compute necessary matrix products
    cx_mat m_1_m_1 = mat[i] * mat[i];

    // compute necessary traces
    double tr_m_1 = trace(mat[i]).real();
    double tr_m_1_m_1 = trace(m_1_m_1).real();

    // off-diagonal update
    if (row_index != column_index) {
      // compute terms D^2 dD^2
      // _______________________________________________________________________________________
      double t_11 = 2 * mat_dim * (m_1_m_1(row_index, row_index).real() + m_1_m_1(column_index, column_index).real());
      double t_21 = 4 * eps[i] * tr_m_1 * (mat[i](row_index, row_index).real() + mat[i](column_index, column_index).real());
      double t_31 = (z * mat[i](column_index, row_index)).real();
      t_31 *= t_31 * 16 * eps[i] * eps[x];
      //________________________________________________________________________________________

      // compute terms D dD D dD
      // _______________________________________________________________________________________

      double t_12 = (mat[i](column_index, row_index) * mat[i](column_index, row_index) * z * z).real();
      t_12 += mat[i](row_index, row_index).real() * mat[i](column_index, column_index).real() * norm(z);
      t_12 *= 4 * mat_dim;

      double t_22 = 4 * eps[i] * tr_m_1 * (mat[i](row_index, row_index).real() + mat[i](column_index, column_index).real());
      double t_32 = (mat[i](column_index, row_index) * z).real();
      t_32 *= t_32 * 16 * eps[i] * eps[x];
      //________________________________________________________________________________________

      // add to total
      temp += 2. * gamma_dim * (norm(z) * (t_11 + t_21 + 4. * tr_m_1_m_1) + t_31);
      temp += cliff * (t_12 + norm(z) * (t_22 + 4. * tr_m_1_m_1) + t_32);
    }

    // diagonal update
    else {
      // compute terms D^2 dD^2
      // _______________________________________________________________________________________
      double t_11 = 2. * mat_dim * m_1_m_1(row_index, row_index).real();
      double t_21 = 4. * eps[x] * m_1_m_1(row_index, row_index).real();
      double t_31 = 4. * eps[i] * tr_m_1 * mat[i](row_index, row_index).real();
      double t_41 = mat[i](row_index, row_index).real();
      t_41 *= t_41 * 4. * eps[i] * eps[x];
      //________________________________________________________________________________________

      // compute terms D dD D dD
      // _______________________________________________________________________________________
      double t_12 = mat[i](row_index, row_index).real();
      t_12 *= t_12 * 2. * mat_dim;
      double t_22 = 4. * eps[x] * m_1_m_1(row_index, row_index).real();
      double t_32 = 4. * eps[i] * tr_m_1 * mat[i](row_index, row_index).real();
      double t_42 = mat[i](row_index, row_index).real();
      t_42 *= t_42 * 4. * eps[i] * eps[x];
      //________________________________________________________________________________________

      // add to total
      temp += 8. * z.real() * z.real() * gamma_dim * (t_11 + t_21 + t_31 + t_41 + 2. * tr_m_1_m_1);
      temp += 4. * z.real() * z.real() * cliff * (t_12 + t_22 + t_32 + t_42 + 2. * tr_m_1_m_1);
    }
  }

  res += 2. * temp;

  // D dD^3 term

  // off-diagonal update
  if (row_index != column_index) {
    temp = 4. * gamma_dim * (mat_dim + 6) * norm(z) * (mat[x](column_index, row_index) * z).real();
    res += 4. * temp;
  }

  // diagonal update
  else {
    double tr_mx = trace(mat[x]).real();
    double rez = 2. * z.real();
    temp = 2. * rez * rez * rez * gamma_dim * (mat[x](row_index, row_index).real() * (mat_dim + 3. * eps[x] + 3.) + eps[x] * tr_mx);
    res += 4. * temp;
  }

  // dD^4 term

  // off-diagonal update
  if (row_index != column_index) {
    temp = gamma_dim * 4. * norm(z) * norm(z) * (mat_dim + 6.);
    res += temp;
  }

  // diagonal update
  else {
    double rez = z.real();
    temp = gamma_dim * 32. * (mat_dim + 3. + 4 * eps[x]) * rez * rez * rez * rez;
    res += temp;
  }

  return res;
}

double Metropolis::runDualAverage(const IDiracOperator& dirac,
                                  const double target) {
  // initial (_i) and final (_f) action2 and action4
  auto* s_i = new double[2];
  auto* s_f = new double[2];
  const auto mat_dim = dirac.getMatrixDimension();
  // calculate length of a sweep in terms of dofs
  const int nsw = dirac.getNumMatrices() * mat_dim * mat_dim - dirac.getNumAntiHermitianMatrices();

  // dual averaging variables
  double stat = 0;
  const double mu = log(10 * m_scale);
  double log_scale_avg = log(m_scale);

  // iter sweeps of metropolis
  for (int i = 0; i < m_num_steps; ++i) {
    for (int j = 0; j < nsw; ++j) {
      constexpr int i_0 = 10;
      constexpr double kappa = 0.75;
      constexpr double shr = 0.05;
      // set action to previous final value,
      // unless it's the first iteration
      if (j) {
        s_i[0] = s_f[0];
        s_i[1] = s_f[1];
      } else {
        s_i[0] = dirac.traceOfDiracSquared();
        s_i[1] = dirac.traceOfDirac4();
      }

      stat += target - runDualAverageCore(dirac, s_i, s_f);

      // perform dual averaging
      const double log_scale = mu - stat * sqrt(i + 1) / (shr * (i + 1 + i_0));
      m_scale = exp(log_scale);
      const double eta = pow(i + 1, -kappa);
      log_scale_avg = eta * log_scale + (1 - eta) * log_scale_avg;
    }
  }

  // set scale on its final dual averaged value
  m_scale = exp(log_scale_avg);

  delete[] s_i;
  delete[] s_f;
  return (stat / (m_num_steps * nsw));
}

double Metropolis::run(const IDiracOperator& dirac) const {
  // initial (_i) and final (_f) action2 and action4
  auto* s_i = new double[2];
  auto* s_f = new double[2];
  const auto mat_dim = dirac.getMatrixDimension();
  // calculate length of a sweep in terms of dofs
  const int nsw = dirac.getNumMatrices() * mat_dim * mat_dim - dirac.getNumAntiHermitianMatrices();

  // return statistic
  double stat = 0;

  // iter sweeps of metropolis
  for (int i = 0; i < m_num_steps; ++i) {
    for (int j = 0; j < nsw; ++j) {
      // set action to previous final value,
      // unless it's the first iteration
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

  delete[] s_i;
  delete[] s_f;

  return (stat / (m_num_steps * nsw));
}

double Metropolis::runDualAverageCore(const IDiracOperator& dirac,
                                      const double* s_i,
                                      double* s_f) const {
  // acceptance probability
  double e;
  const auto num_matrices = dirac.getNumMatrices();
  const auto mat_dim = dirac.getMatrixDimension();

  // metropolis
  const int x = (int)(num_matrices * m_rng->getUniform());
  const int row_index = (int)(mat_dim * m_rng->getUniform());
  const int column_index = (int)(mat_dim * m_rng->getUniform());

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

  const double delta_2 = delta2(dirac, x, row_index, column_index, z);
  const double delta_4 = delta4(dirac, x, row_index, column_index, z);
  const double action_delta = delta24(dirac, x, row_index, column_index, z);

  auto& mat = dirac.getMatrices();
  // metropolis test
  if (action_delta < 0) {
    // update matrix element
    if (row_index != column_index) {
      mat[x](row_index, column_index) += z;
      mat[x](column_index, row_index) += conj(z);
    } else {
      mat[x](row_index, row_index) += 2. * z;
    }

    // update action
    s_f[0] = s_i[0] + delta_2;
    s_f[1] = s_i[1] + delta_4;

    // move accepted
    e = 1;
  } else {
    e = exp(-action_delta);
    const double p = m_rng->getUniform();

    if (e > p) {
      // update matrix element
      if (row_index != column_index) {
        mat[x](row_index, column_index) += z;
        mat[x](column_index, row_index) += conj(z);
      } else {
        mat[x](row_index, row_index) += 2. * z;
      }

      // update action
      s_f[0] = s_i[0] + delta_2;
      s_f[1] = s_i[1] + delta_4;
    } else {
      s_f[0] = s_i[0];
      s_f[1] = s_i[1];
    }
  }

  return e;
}

double Metropolis::runCore(const IDiracOperator& dirac,
                           const double* s_i,
                           double* s_f) const {
  // acceptance probability
  double ret = 0;

  const auto num_matrices = dirac.getNumMatrices();
  const auto mat_dim = dirac.getMatrixDimension();

  // metropolis
  const int x = (int)(num_matrices * m_rng->getUniform());
  const int row_index = (int)(mat_dim * m_rng->getUniform());
  const int column_index = (int)(mat_dim * m_rng->getUniform());
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

  const double delta_2 = delta2(dirac, x, row_index, column_index, z);
  const double delta_4 = delta4(dirac, x, row_index, column_index, z);
  const double action_delta = m_action->getG2() * delta_2 + delta_4;

  auto& mat = dirac.getMatrices();
  // metropolis test
  if (action_delta < 0) {
    // update matrix element
    if (row_index != column_index) {
      mat[x](row_index, column_index) += z;
      mat[x](column_index, row_index) += conj(z);
    } else {
      mat[x](row_index, row_index) += 2. * z;
    }

    // update action
    s_f[0] = s_i[0] + delta_2;
    s_f[1] = s_i[1] + delta_4;

    // move accepted
    ret = 1;
  } else {
    const double e = exp(-action_delta);
    const double p = m_rng->getUniform();

    if (e > p) {
      // update matrix element
      if (row_index != column_index) {
        mat[x](row_index, column_index) += z;
        mat[x](column_index, row_index) += conj(z);
      } else {
        mat[x](row_index, row_index) += 2. * z;
      }

      // update action
      s_f[0] = s_i[0] + delta_2;
      s_f[1] = s_i[1] + delta_4;

      // move accepted
      ret = 1;
    } else {
      s_f[0] = s_i[0];
      s_f[1] = s_i[1];
    }
  }

  return ret;
}
