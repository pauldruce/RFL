//
// Created by Paul Druce on 17/11/2023.
//

#include "IDiracOperatorDerivatives.hpp"

using namespace arma;

static cx_mat computeB4(const IDiracOperator& dirac, const int& k,
                        const int& i_2,
                        const int& i_3,
                        const int& i_4,
                        const double& cliff,
                        const bool& neg) {
  auto& matrices = dirac.getMatrices();
  auto& epsilons = dirac.getEpsilons();
  const auto dim = dirac.getMatrixDimension();

  // base matrix products
  cx_mat m_2_m_3 = matrices[i_2] * matrices[i_3];
  cx_mat m_2_m_4 = matrices[i_2] * matrices[i_4];
  cx_mat m_3_m_4 = matrices[i_3] * matrices[i_4];
  cx_mat m_2_m_3_m_4 = m_2_m_3 * matrices[i_4];

  // return value
  cx_mat res(dim, dim, fill::eye);

  if (neg) {

    // traces
    double tr_234 = trace(m_2_m_3_m_4).imag();
    double tr_2 = trace(matrices[i_2]).real();
    double tr_3 = trace(matrices[i_3]).real();
    double tr_4 = trace(matrices[i_4]).real();

    // compute sum
    res *= -2 * epsilons[k] * tr_234;
    res += cx_double(0., dim) * (m_2_m_3_m_4 - m_2_m_3_m_4.t());
    res += cx_double(0., epsilons[i_2] * tr_2) * (m_3_m_4 - m_3_m_4.t());
    res += cx_double(0., epsilons[i_3] * tr_3) * (m_2_m_4 - m_2_m_4.t());
    res += cx_double(0., epsilons[i_4] * tr_4) * (m_2_m_3 - m_2_m_3.t());
  } else {
    // traces
    double tr_234 = trace(m_2_m_3_m_4).real();
    double tr_23 = trace(m_2_m_3).real();
    double tr_24 = trace(m_2_m_4).real();
    double tr_34 = trace(m_3_m_4).real();
    double tr_2 = trace(matrices[i_2]).real();
    double tr_3 = trace(matrices[i_3]).real();
    double tr_4 = trace(matrices[i_4]).real();

    // compute sum
    res *= 2 * epsilons[k] * tr_234;
    res += dim * (m_2_m_3_m_4 + m_2_m_3_m_4.t());
    res += epsilons[i_2] * tr_2 * (m_3_m_4 + m_3_m_4.t());
    res += epsilons[i_3] * tr_3 * (m_2_m_4 + m_2_m_4.t());
    res += epsilons[i_4] * tr_4 * (m_2_m_3 + m_2_m_3.t());
    res += 2 * epsilons[k] * epsilons[i_2] * tr_34 * matrices[i_2];
    res += 2 * epsilons[k] * epsilons[i_3] * tr_24 * matrices[i_3];
    res += 2 * epsilons[k] * epsilons[i_4] * tr_23 * matrices[i_4];
  }

  return cliff * res;
}

static cx_mat computeB2(const IDiracOperator& dirac, const int& k, const int& i) {
  const auto num_matrices = dirac.getNumMatrices();
  const auto dim = dirac.getMatrixDimension();
  const auto gamma_dim = dirac.getGammaDimension();

  auto& omega_table_4 = dirac.getOmegaTable4();
  auto& matrices = dirac.getMatrices();
  auto& epsilons = dirac.getEpsilons();

  // clifford product
  double cliff = omega_table_4[i + num_matrices * (k + num_matrices * (i + num_matrices * k))].real();

  // base matrix products
  cx_mat mi_mk = matrices[i] * matrices[k];
  cx_mat mi_mi = matrices[i] * matrices[i];
  cx_mat mi_mi_mk = matrices[i] * mi_mk;
  cx_mat mi_mk_mi = mi_mk * matrices[i];

  // traces
  double triki = trace(mi_mk_mi).real();
  double trik = trace(mi_mk).real();
  double trii = trace(mi_mi).real();
  double tri = trace(matrices[i]).real();
  double trk = trace(matrices[k]).real();

  // return value
  cx_mat res(dim, dim, fill::eye);

  if (cliff < 0) {
    // compute sum
    res *= epsilons[k] * triki;
    res += dim * (mi_mi_mk + mi_mi_mk.t() - mi_mk_mi);
    res += epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 2 * epsilons[k] * epsilons[i] * trik * matrices[i];
    res += epsilons[k] * trk * mi_mi;
    res += trii * matrices[k];
  } else {
    // compute sum
    res *= 3 * epsilons[k] * triki;
    res += dim * (mi_mi_mk + mi_mi_mk.t() + mi_mk_mi);
    res += 3 * epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 6 * epsilons[k] * epsilons[i] * trik * matrices[i];
    res += 3 * epsilons[k] * trk * mi_mi;
    res += 3 * trii * matrices[k];
  }

  return 2 * gamma_dim * res;
}

static cx_mat computeB(const IDiracOperator& dirac, const int& k) {
  const auto m_dim = dirac.getMatrixDimension();
  const auto m_gamma_dim = dirac.getGammaDimension();
  const auto& matrices = dirac.getMatrices();
  const auto& epsilons = dirac.getEpsilons();

  // base matrix products
  const cx_mat m_2 = matrices[k] * matrices[k];
  const cx_mat m_3 = matrices[k] * m_2;

  // traces
  const double tr_3 = trace(m_3).real();
  const double tr_2 = trace(m_2).real();
  const double tr_1 = trace(matrices[k]).real();

  cx_mat res(m_dim, m_dim, fill::eye);
  res *= epsilons[k] * tr_3;
  res += m_dim * m_3;
  res += 3 * tr_2 * matrices[k];
  res += 3 * epsilons[k] * tr_1 * m_2;

  return 2 * m_gamma_dim * res;
}

cx_mat derDirac24(const IDiracOperator& dirac, const int& k, const bool& herm, const double g_2) {
  return g_2 * derDirac2(dirac, k) + derDirac4(dirac, k, herm);
}

arma::cx_mat derDirac2(const IDiracOperator& dirac, const int& k) {
  const auto dim = dirac.getMatrixDimension();
  const auto gamma_dim = dirac.getGammaDimension();
  const auto& epsilons = dirac.getEpsilons();
  const auto& matrices = dirac.getMatrices();
  arma::cx_mat res(dim, dim, fill::eye);
  res *= epsilons[k] * trace(matrices[k]).real();
  res += dim * matrices[k];
  return 4 * gamma_dim * res;
}

cx_mat derDirac4(const IDiracOperator& dirac, const int& k, const bool& herm) {
  const auto dim = dirac.getMatrixDimension();
  const auto num_matrices = dirac.getNumMatrices();
  const auto& epsilons = dirac.getEpsilons();
  const auto& omega_table_4 = dirac.getOmegaTable4();

  cx_mat res(dim, dim, fill::zeros);

  // four distinct indices
  for (int i_1 = 0; i_1 < num_matrices; ++i_1) {
    if (i_1 != k) {
      for (int i_2 = i_1 + 1; i_2 < num_matrices; ++i_2) {
        if (i_2 != k) {
          for (int i_3 = i_2 + 1; i_3 < num_matrices; ++i_3) {
            if (i_3 != k) {
              // epsilon factor

              if (double e = epsilons[k] * epsilons[i_1] * epsilons[i_2] * epsilons[i_3]; e < 0) {
                // clifford product
                double cliff_1 = omega_table_4[i_3 + num_matrices * (i_2 + num_matrices * (i_1 + num_matrices * k))].imag();
                double cliff_2 = omega_table_4[i_2 + num_matrices * (i_3 + num_matrices * (i_1 + num_matrices * k))].imag();
                double cliff_3 = omega_table_4[i_3 + num_matrices * (i_1 + num_matrices * (i_2 + num_matrices * k))].imag();

                if (fabs(cliff_1) > 1e-10) {
                  res += computeB4(dirac, k, i_1, i_2, i_3, cliff_1, true);
                  res += computeB4(dirac, k, i_1, i_3, i_2, cliff_2, true);
                  res += computeB4(dirac, k, i_2, i_1, i_3, cliff_3, true);
                }
              } else {
                // clifford product
                double cliff_1 = omega_table_4[i_3 + num_matrices * (i_2 + num_matrices * (i_1 + num_matrices * k))].real();
                double cliff_2 = omega_table_4[i_2 + num_matrices * (i_3 + num_matrices * (i_1 + num_matrices * k))].real();
                double cliff_3 = omega_table_4[i_3 + num_matrices * (i_1 + num_matrices * (i_2 + num_matrices * k))].real();

                if (fabs(cliff_1) > 1e-10) {
                  res += computeB4(dirac, k, i_1, i_2, i_3, cliff_1, false);
                  res += computeB4(dirac, k, i_1, i_3, i_2, cliff_2, false);
                  res += computeB4(dirac, k, i_2, i_1, i_3, cliff_3, false);
                }
              }
            }
          }
        }
      }
    }
  }
  res = res + res.t();

  // two distinct pairs of equal indices
  for (int i = 0; i < num_matrices; ++i) {
    if (i != k) {
      res += computeB2(dirac, k, i);
    }
  }

  // all indices equal
  res += computeB(dirac, k);

  if (herm) {
    return 2 * (res + res.t());
  } else {
    return 4 * res;
  }
}
