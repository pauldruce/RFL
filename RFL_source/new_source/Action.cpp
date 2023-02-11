//
// Created by Paul Druce on 12/11/2022.
//
#include "Action.hpp"

using namespace std;
using namespace arma;

Action::Action(double g_2) : m_g_2(g_2), m_g_4(1.0) { }

Action::Action(double g_2, double g_4) : m_g_2(g_2), m_g_4(g_4) { }

double Action::calculateSFromDirac(const DiracOperator& dirac) const {
  cx_mat dirac_mat = dirac.getDiracMatrix();
  cx_mat dirac_squared = dirac_mat * dirac_mat;
  double trace_dirac_squared = trace(dirac_squared).real();
  double trace_dirac_4 = trace(dirac_squared * dirac_squared).real();
  return m_g_2 * trace_dirac_squared + m_g_4 * trace_dirac_4;
}

double Action::calculateS(const DiracOperator& dirac) const {
  return m_g_2 * dirac2(dirac) + m_g_4 * dirac4(dirac);
}

double Action::dirac2(const DiracOperator& dirac) const {
  auto* mat = dirac.getMatrices();
  auto* eps = dirac.getEpsilons();
  double res = 0.;
  for (int i = 0; i < dirac.getNumMatrices(); ++i) {
    double trace_squared = trace(mat[i] * mat[i]).real();
    double trace_mat = trace(mat[i]).real();

    res += (dirac.getMatrixDimension() * trace_squared + eps[i] * trace_mat * trace_mat);
  }

  return 2. * dirac.getGammaDimension() * res;
}

double Action::dirac4(const DiracOperator& dirac) const {
  double res = 0.;

  auto num_matrices = dirac.getNumMatrices();

  // four distinct indices
  for (int i = 0; i < num_matrices; ++i) {
    for (int j = i + 1; j < num_matrices; ++j) {
      for (int k = j + 1; k < num_matrices; ++k) {
        for (int l = k + 1; l < num_matrices; ++l)
          res += 8 * (computeA4(dirac, i, j, k, l) + computeA4(dirac, i, j, l, k) + computeA4(dirac, i, k, j, l));
      }
    }
  }

  // two distinct pairs of equal indices
  for (int i = 0; i < num_matrices; ++i) {
    for (int j = i + 1; j < num_matrices; ++j)
      res += 2 * computeA2(dirac, i, j);
  }

  // all indices equal
  for (int i = 0; i < num_matrices; ++i)
    res += computeA(dirac, i);

  return res;
}

double Action::computeA4(const DiracOperator& dirac, const int& i_1, const int& i_2, const int& i_3, const int& i_4) const {
  // epsilon factor
  auto* eps = dirac.getEpsilons();
  auto* omega_table_4 = dirac.getOmegaTable4();
  auto* mat = dirac.getMatrices();
  auto num_matrices = dirac.getNumMatrices();
  auto mat_dim = dirac.getMatrixDimension();

  const int e = eps[i_1] * eps[i_2] * eps[i_3] * eps[i_4];

  // if e=-1, then [1+*e] becomes 2i*imag
  // and the clifford part is guaranteed to
  // be pure imaginary
  if (e < 0) {
    // clifford product
    const double cliff = omega_table_4[i_4 + num_matrices * (i_3 + num_matrices * (i_2 + num_matrices * i_1))].imag();

    if (fabs(cliff) > 1e-10) {
      // base matrix products
      const cx_mat m_1_m_2 = mat[i_1] * mat[i_2];
      const cx_mat m_1_m_3 = mat[i_1] * mat[i_3];
      const cx_mat m_1_m_4 = mat[i_1] * mat[i_4];
      const cx_mat m_2_m_3 = mat[i_2] * mat[i_3];
      const cx_mat m_2_m_4 = mat[i_2] * mat[i_4];
      const cx_mat m_3_m_4 = mat[i_3] * mat[i_4];

      // traces
      double tr_1234 = trace(m_1_m_2 * m_3_m_4).imag();
      double tr_234 = trace(m_2_m_3 * mat[i_4]).imag();
      double tr_134 = trace(m_1_m_3 * mat[i_4]).imag();
      double tr_124 = trace(m_1_m_2 * mat[i_4]).imag();
      double tr_123 = trace(m_1_m_2 * mat[i_3]).imag();
      double tr_1 = trace(mat[i_1]).real();
      double tr_2 = trace(mat[i_2]).real();
      double tr_3 = trace(mat[i_3]).real();
      double tr_4 = trace(mat[i_4]).real();

      // compute sum
      double res = mat_dim * tr_1234;
      res += eps[i_1] * tr_1 * tr_234;
      res += eps[i_2] * tr_2 * tr_134;
      res += eps[i_3] * tr_3 * tr_124;
      res += eps[i_4] * tr_4 * tr_123;

      return -2 * cliff * res;
      // NOTE: this minus here comes from the 'i' in cliff
      // and the 'i' coming from 2i*imag
    } else {
      return 0.;
    }
  } else {
    // clifford product
    const double cliff = omega_table_4[i_4 + num_matrices * (i_3 + num_matrices * (i_2 + num_matrices * i_1))].real();

    if (fabs(cliff) > 1e-10) {
      // base matrix products
      cx_mat m_1_m_2 = mat[i_1] * mat[i_2];
      cx_mat m_1_m_3 = mat[i_1] * mat[i_3];
      cx_mat m_1_m_4 = mat[i_1] * mat[i_4];
      cx_mat m_2_m_3 = mat[i_2] * mat[i_3];
      cx_mat m_2_m_4 = mat[i_2] * mat[i_4];
      cx_mat m_3_m_4 = mat[i_3] * mat[i_4];

      // traces
      double tr_1234 = trace(m_1_m_2 * m_3_m_4).real();
      double tr_234 = trace(m_2_m_3 * mat[i_4]).real();
      double tr_134 = trace(m_1_m_3 * mat[i_4]).real();
      double tr_124 = trace(m_1_m_2 * mat[i_4]).real();
      double tr_123 = trace(m_1_m_2 * mat[i_3]).real();
      double tr_12 = trace(m_1_m_2).real();
      double tr_34 = trace(m_3_m_4).real();
      double tr_13 = trace(m_1_m_3).real();
      double tr_24 = trace(m_2_m_4).real();
      double tr_14 = trace(m_1_m_4).real();
      double tr_23 = trace(m_2_m_3).real();
      double tr_1 = trace(mat[i_1]).real();
      double tr_2 = trace(mat[i_2]).real();
      double tr_3 = trace(mat[i_3]).real();
      double tr_4 = trace(mat[i_4]).real();

      double res = mat_dim * tr_1234;
      res += eps[i_1] * tr_1 * tr_234;
      res += eps[i_2] * tr_2 * tr_134;
      res += eps[i_3] * tr_3 * tr_124;
      res += eps[i_4] * tr_4 * tr_123;
      res += eps[i_1] * eps[i_2] * tr_12 * tr_34;
      res += eps[i_1] * eps[i_3] * tr_13 * tr_24;
      res += eps[i_1] * eps[i_4] * tr_14 * tr_23;

      //	  cliff = omega_table_4[i4 + D.nHL * (i3 + D.nHL * (i2 + D.nHL * i1))].real();

      return 2 * cliff * res;
    } else {
      return 0.;
    }
  }
}

double Action::computeA2(const DiracOperator& dirac, const int& i_1, const int& i_2) const {
  auto* omega_table_4 = dirac.getOmegaTable4();
  auto* mat = dirac.getMatrices();
  auto* eps = dirac.getEpsilons();
  auto num_matrices = dirac.getNumMatrices();
  auto mat_dim = dirac.getMatrixDimension();
  auto gamma_dim = dirac.getGammaDimension();

  // clifford product
  double cliff = omega_table_4[i_2 + num_matrices * (i_1 + num_matrices * (i_2 + num_matrices * i_1))].real();

  // base matrix products
  cx_mat m_1_m_1 = mat[i_1] * mat[i_1];
  cx_mat m_2_m_2 = mat[i_2] * mat[i_2];
  cx_mat m_1_m_2 = mat[i_1] * mat[i_2];

  // traces
  double tr_1122 = trace(m_1_m_1 * m_2_m_2).real();
  double tr_1212 = trace(m_1_m_2 * m_1_m_2).real();
  double tr_122 = trace(m_1_m_2 * mat[i_2]).real();
  double tr_112 = trace(m_1_m_1 * mat[i_2]).real();
  double tr_11 = trace(m_1_m_1).real();
  double tr_22 = trace(m_2_m_2).real();
  double tr_12 = trace(m_1_m_2).real();
  double tr_1 = trace(mat[i_1]).real();
  double tr_2 = trace(mat[i_2]).real();

  if (cliff < 0) {
    // compute sum
    double res = mat_dim * (2 * tr_1122 - tr_1212);
    res += 2 * eps[i_1] * tr_1 * tr_122;
    res += 2 * eps[i_2] * tr_2 * tr_112;
    res += tr_11 * tr_22;
    res += 2 * eps[i_1] * eps[i_2] * tr_12 * tr_12;

    return 2 * gamma_dim * res;
  } else {
    // compute sum
    double res = mat_dim * (2 * tr_1122 + tr_1212);
    res += 6 * eps[i_1] * tr_1 * tr_122;
    res += 6 * eps[i_2] * tr_2 * tr_112;
    res += 3 * tr_11 * tr_22;
    res += 6 * eps[i_1] * eps[i_2] * tr_12 * tr_12;

    return 2 * gamma_dim * res;
  }
}

double Action::computeA(const DiracOperator& dirac, const int& i) const {
  auto* mat = dirac.getMatrices();
  auto* eps = dirac.getEpsilons();
  auto mat_dim = dirac.getMatrixDimension();
  auto gamma_dim = dirac.getGammaDimension();

  // base matrix products
  cx_mat m_2 = mat[i] * mat[i];
  cx_mat m_3 = m_2 * mat[i];

  // traces
  double tr_1 = trace(mat[i]).real();
  double tr_2 = trace(m_2).real();
  double tr_3 = trace(m_3).real();
  double tr_4 = trace(m_3 * mat[i]).real();

  double res = mat_dim * tr_4;
  res += 4 * eps[i] * tr_1 * tr_3;
  res += 3 * tr_2 * tr_2;

  return 2 * gamma_dim * res;
}

void Action::setParams(double g_2, double g_4) {
  this->m_g_2 = g_2;
  this->m_g_4 = g_4;
}
void Action::setG4(double value) { this->m_g_4 = value; }
void Action::setG2(double value) { this->m_g_2 = value; }
