//
// Implementation of BarrettGlaserAction.
//

#include "BarrettGlaserAction.hpp"
#include <cmath>

using namespace std;
using namespace arma;

BarrettGlaserAction::BarrettGlaserAction(const double g_2) : m_g_2(g_2), m_g_4(1.0) {}

BarrettGlaserAction::BarrettGlaserAction(const double g_2, const double g_4) : m_g_2(g_2), m_g_4(g_4) {}

void BarrettGlaserAction::setG2(const double value) { this->m_g_2 = value; }

void BarrettGlaserAction::setG4(const double value) { this->m_g_4 = value; }

void BarrettGlaserAction::setParams(const double g_2, const double g_4) {
  this->m_g_2 = g_2;
  this->m_g_4 = g_4;
}

// Evaluates S(D) directly from powers of the assembled matrix D:
// D^2 = D * D, Tr(D^2) = Re(tr(D^2)), and Tr(D^4) = Re(tr(D^2 * D^2)).
// S(D) = g_2 Tr(D^2) + g_4 Tr(D^4).
double BarrettGlaserAction::calculateSFromDirac(const IDiracOperator& dirac) const {
  const cx_mat dirac_mat = dirac.getDiracMatrix();
  const cx_mat dirac_squared = dirac_mat * dirac_mat;
  const double trace_dirac_squared = trace(dirac_squared).real();
  const double trace_dirac_4 = trace(dirac_squared * dirac_squared).real();
  return m_g_2 * trace_dirac_squared + m_g_4 * trace_dirac_4;
}

// Evaluates S(D) using component matrix traces from IDiracOperator:
// S(D) = g_2 Tr(D^2) + g_4 Tr(D^4).
// Avoids assembling the full (d_gamma * N) x (d_gamma * N) matrix D.
double BarrettGlaserAction::calculateS(const IDiracOperator& dirac) const {
  return m_g_2 * dirac.traceOfDiracSquared() + m_g_4 * dirac.traceOfDirac4();
}

// Evaluates the exact analytic variation Delta Tr(D^2) under rank-1 Hermitian move:
// delta_M = z * E_{ij} + conj(z) * E_{ji}.
// - Off-diagonal (i != j): Delta_2 = 4 * d_gamma * N * (2 * Re(z * M(j, i)) + |z|^2).
// - Diagonal (i == j): delta_M_ii = 2 * Re(z), yielding:
//   Delta_2 = 8 * d_gamma * Re(z) * [ N (M(i, i) + Re(z)) + eps * (Tr(M) + Re(z)) ].
double BarrettGlaserAction::delta2(
    const IDiracOperator& dirac,
    const int x,
    const int row_index,
    const int column_index,
    const cx_double& z) {
  const auto& mat = dirac.getMatrices();
  const auto& eps = dirac.getEpsilons();
  const auto mat_dim = dirac.getMatrixDimension();
  const auto gamma_dim = dirac.getGammaDimension();

  if (row_index != column_index) {
    // Off-diagonal variation: delta_M = z * E_{ij} + conj(z) * E_{ji}.
    // Hermiticity requires perturbing both (i, j) and (j, i) simultaneously.
    return 4. * gamma_dim * mat_dim * (2. * (z * mat[x](column_index, row_index)).real() + norm(z));
  } else {
    // Diagonal variation: delta_M = 2 * Re(z) * E_{ii}.
    // Preserves Hermiticity by requiring Im(z) = 0; the symmetric perturbation collapses to 2 * Re(z).
    const double tr_m = trace(mat[x]).real();
    return 8. * gamma_dim * z.real() * (mat_dim * (mat[x](row_index, row_index).real() + z.real()) + eps[x] * (tr_m + z.real()));
  }
}

// Evaluates the exact analytic variation Delta Tr(D^4) under delta D = gamma^x (x) delta_M.
// Binomial expansion yields terms of order 1 to 4 in perturbation:
// Delta_4 = 4 Tr(D^3 dD) + 2 Tr(D^2 dD^2) + Tr(D dD D dD) + 4 Tr(D dD^3) + Tr(dD^4).
// Evaluated in O(N^2) operations using precomputed Clifford trace tensor Omega.
double BarrettGlaserAction::delta4(
    const IDiracOperator& dirac,
    const int x,
    const int row_index,
    const int column_index,
    const cx_double& z) {
  double res = 0.;

  const auto& omega_table_4 = dirac.getOmegaTable4();
  const auto& mat = dirac.getMatrices();
  const auto& eps = dirac.getEpsilons();
  const auto mat_dim = dirac.getMatrixDimension();
  const auto gamma_dim = dirac.getGammaDimension();
  const auto num_matrices = dirac.getNumMatrices();

  // Linear variation term: 4 Tr(D^3 dD).
  // Symmetries of the Clifford trace tensor Omega allow index folding:
  // folding the sum over i_1 <= i_3 introduces a factor of 2 for distinct indices (i_1 != i_3).
  for (int i_3 = 0; i_3 < num_matrices; ++i_3) {
    for (int i_2 = 0; i_2 < num_matrices; ++i_2) {
      for (int i_1 = 0; i_1 <= i_3; ++i_1) {
        cx_double cliff = omega_table_4[x + num_matrices * (i_3 + num_matrices * (i_2 + num_matrices * i_1))];

        if (std::abs(cliff.real()) > 1e-10 || std::abs(cliff.imag()) > 1e-10) {
          // Compute necessary matrix products and traces.
          cx_mat m_1_m_2 = mat[i_1] * mat[i_2];
          cx_mat m_2_m_3 = mat[i_2] * mat[i_3];
          cx_mat m_1_m_3 = mat[i_1] * mat[i_3];
          cx_mat m_1_m_2_m_3 = mat[i_1] * m_2_m_3;

          double tr_m_1 = trace(mat[i_1]).real();
          double tr_m_2 = trace(mat[i_2]).real();
          double tr_m_3 = trace(mat[i_3]).real();
          double tr_m_1_m_2 = trace(m_1_m_2).real();
          double tr_m_2_m_3 = trace(m_2_m_3).real();
          double tr_m_1_m_3 = trace(m_1_m_3).real();
          cx_double tr_m_1_m_2_m_3 = trace(m_1_m_2_m_3);

          if (row_index != column_index) {
            // Off-diagonal variation: delta_M = z * E_{ij} + conj(z) * E_{ji}.
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

            // Fold symmetric indices: distinct pairs (i_1 != i_3) carry multiplicity 2.
            if (i_1 != i_3) {
              res += 2. * (cliff * t_1).real();
            } else {
              res += (cliff * t_1).real();
            }
          } else {
            // Diagonal variation: delta_M = 2 * Re(z) * E_{ii}.
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

            // Fold symmetric indices: distinct pairs (i_1 != i_3) carry multiplicity 2.
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

  // Quadratic variation terms: 2 Tr(D^2 dD^2) + Tr(D dD D dD).
  // Cross terms between Dirac operator components and perturbation matrices.
  double temp = 0;
  for (int i = 0; i < num_matrices; ++i) {
    double cliff = omega_table_4[x + num_matrices * (i + num_matrices * (x + num_matrices * i))].real();

    cx_mat m_1_m_1 = mat[i] * mat[i];

    double tr_m_1 = trace(mat[i]).real();
    double tr_m_1_m_1 = trace(m_1_m_1).real();

    if (row_index != column_index) {
      double t_11 = 2.0 * mat_dim * (m_1_m_1(row_index, row_index).real() + m_1_m_1(column_index, column_index).real());
      double t_21 = 4.0 * eps[i] * tr_m_1 * (mat[i](row_index, row_index).real() + mat[i](column_index, column_index).real());
      double t_31 = (z * mat[i](column_index, row_index)).real();
      t_31 *= t_31 * 16.0 * eps[i] * eps[x];

      double t_12 = (mat[i](column_index, row_index) * mat[i](column_index, row_index) * z * z).real();
      t_12 += mat[i](row_index, row_index).real() * mat[i](column_index, column_index).real() * norm(z);
      t_12 *= 4.0 * mat_dim;

      double t_22 = 4.0 * eps[i] * tr_m_1 * (mat[i](row_index, row_index).real() + mat[i](column_index, column_index).real());
      double t_32 = (mat[i](column_index, row_index) * z).real();
      t_32 *= t_32 * 16.0 * eps[i] * eps[x];

      temp += 2.0 * gamma_dim * (norm(z) * (t_11 + t_21 + 4.0 * tr_m_1_m_1) + t_31);
      temp += cliff * (t_12 + norm(z) * (t_22 + 4.0 * tr_m_1_m_1) + t_32);
    } else {
      double t_11 = 2.0 * mat_dim * m_1_m_1(row_index, row_index).real();
      double t_21 = 4.0 * eps[x] * m_1_m_1(row_index, row_index).real();
      double t_31 = 4.0 * eps[i] * tr_m_1 * mat[i](row_index, row_index).real();
      double t_41 = mat[i](row_index, row_index).real();
      t_41 *= t_41 * 4.0 * eps[i] * eps[x];

      double t_12 = mat[i](row_index, row_index).real();
      t_12 *= t_12 * 2.0 * mat_dim;
      double t_22 = 4.0 * eps[x] * m_1_m_1(row_index, row_index).real();
      double t_32 = 4.0 * eps[i] * tr_m_1 * mat[i](row_index, row_index).real();
      double t_42 = mat[i](row_index, row_index).real();
      t_42 *= t_42 * 4.0 * eps[i] * eps[x];

      temp += 8.0 * z.real() * z.real() * gamma_dim * (t_11 + t_21 + t_31 + t_41 + 2.0 * tr_m_1_m_1);
      temp += 4.0 * z.real() * z.real() * cliff * (t_12 + t_22 + t_32 + t_42 + 2.0 * tr_m_1_m_1);
    }
  }

  res += 2. * temp;

  // Cubic variation term: 4 Tr(D dD^3).
  if (row_index != column_index) {
    temp = 4. * gamma_dim * (mat_dim + 6) * norm(z) * (mat[x](column_index, row_index) * z).real();
    res += 4. * temp;
  } else {
    double tr_mx = trace(mat[x]).real();
    double rez = 2. * z.real();
    temp = 2. * rez * rez * rez * gamma_dim * (mat[x](row_index, row_index).real() * (mat_dim + 3. * eps[x] + 3.) + eps[x] * tr_mx);
    res += 4. * temp;
  }

  // Quartic variation term: Tr(dD^4).
  if (row_index != column_index) {
    temp = gamma_dim * 4. * norm(z) * norm(z) * (mat_dim + 6.);
    res += temp;
  } else {
    double rez = z.real();
    temp = gamma_dim * 32. * (mat_dim + 3. + 4 * eps[x]) * rez * rez * rez * rez;
    res += temp;
  }

  return res;
}

// Combines quadratic and quartic trace variations linearly:
// Delta S = g_2 * Delta_2 + g_4 * Delta_4.
double BarrettGlaserAction::delta24(
    const IDiracOperator& dirac,
    const int x,
    const int row_index,
    const int column_index,
    const cx_double& z,
    const double g_2,
    const double g_4) {
  return g_2 * delta2(dirac, x, row_index, column_index, z) + g_4 * delta4(dirac, x, row_index, column_index, z);
}

double BarrettGlaserAction::delta24(
    const IDiracOperator& dirac,
    const int x,
    const int row_index,
    const int column_index,
    const cx_double& z) const {
  return delta24(dirac, x, row_index, column_index, z, m_g_2, m_g_4);
}

// ----------------------------------------------------------------------------
// Analytical Gradient Helpers
// ----------------------------------------------------------------------------

namespace {

// Evaluates the 4-distinct-index contraction contribution to grad(Tr(D^4)):
// B4 accounts for configurations where indices (k, i_2, i_3, i_4) are all distinct.
// When neg is true, the Clifford trace is purely imaginary and contracts anti-symmetrically.
cx_mat computeB4(
    const IDiracOperator& dirac,
    const int k,
    const int i_2,
    const int i_3,
    const int i_4,
    const double cliff,
    const bool neg) {
  const auto& matrices = dirac.getMatrices();
  const auto& epsilons = dirac.getEpsilons();
  const auto dim = dirac.getMatrixDimension();

  cx_mat m_2_m_3 = matrices[i_2] * matrices[i_3];
  cx_mat m_2_m_4 = matrices[i_2] * matrices[i_4];
  cx_mat m_3_m_4 = matrices[i_3] * matrices[i_4];
  cx_mat m_2_m_3_m_4 = m_2_m_3 * matrices[i_4];

  cx_mat res(dim, dim, fill::eye);

  if (neg) {
    double tr_234 = trace(m_2_m_3_m_4).imag();
    double tr_2 = trace(matrices[i_2]).real();
    double tr_3 = trace(matrices[i_3]).real();
    double tr_4 = trace(matrices[i_4]).real();

    res *= -2 * epsilons[k] * tr_234;
    res += cx_double(0., dim) * (m_2_m_3_m_4 - m_2_m_3_m_4.t());
    res += cx_double(0., epsilons[i_2] * tr_2) * (m_3_m_4 - m_3_m_4.t());
    res += cx_double(0., epsilons[i_3] * tr_3) * (m_2_m_4 - m_2_m_4.t());
    res += cx_double(0., epsilons[i_4] * tr_4) * (m_2_m_3 - m_2_m_3.t());
  } else {
    double tr_234 = trace(m_2_m_3_m_4).real();
    double tr_23 = trace(m_2_m_3).real();
    double tr_24 = trace(m_2_m_4).real();
    double tr_34 = trace(m_3_m_4).real();
    double tr_2 = trace(matrices[i_2]).real();
    double tr_3 = trace(matrices[i_3]).real();
    double tr_4 = trace(matrices[i_4]).real();

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

// Evaluates the 2-distinct-index contraction contribution to grad(Tr(D^4)):
// B2 accounts for index coincidences where two matrices share index i and two share index k.
cx_mat computeB2(const IDiracOperator& dirac, const int k, const int i) {
  const auto num_matrices = dirac.getNumMatrices();
  const auto dim = dirac.getMatrixDimension();
  const auto gamma_dim = dirac.getGammaDimension();

  const auto& omega_table_4 = dirac.getOmegaTable4();
  const auto& matrices = dirac.getMatrices();
  const auto& epsilons = dirac.getEpsilons();

  double cliff = omega_table_4[i + num_matrices * (k + num_matrices * (i + num_matrices * k))].real();

  cx_mat mi_mk = matrices[i] * matrices[k];
  cx_mat mi_mi = matrices[i] * matrices[i];
  cx_mat mi_mi_mk = matrices[i] * mi_mk;
  cx_mat mi_mk_mi = mi_mk * matrices[i];

  double triki = trace(mi_mk_mi).real();
  double trik = trace(mi_mk).real();
  double trii = trace(mi_mi).real();
  double tri = trace(matrices[i]).real();
  double trk = trace(matrices[k]).real();

  cx_mat res(dim, dim, fill::eye);

  if (cliff < 0) {
    res *= epsilons[k] * triki;
    res += dim * (mi_mi_mk + mi_mi_mk.t() - mi_mk_mi);
    res += epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 2 * epsilons[k] * epsilons[i] * trik * matrices[i];
    res += epsilons[k] * trk * mi_mi;
    res += trii * matrices[k];
  } else {
    res *= 3 * epsilons[k] * triki;
    res += dim * (mi_mi_mk + mi_mi_mk.t() + mi_mk_mi);
    res += 3 * epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 6 * epsilons[k] * epsilons[i] * trik * matrices[i];
    res += 3 * epsilons[k] * trk * mi_mi;
    res += 3 * trii * matrices[k];
  }

  return 2 * gamma_dim * res;
}

// Evaluates the single-index self-coupling contribution to grad(Tr(D^4)):
// B accounts for the quartic trace term Tr(M_k^4) when all four indices equal k.
cx_mat computeB(const IDiracOperator& dirac, const int k) {
  const auto m_dim = dirac.getMatrixDimension();
  const auto m_gamma_dim = dirac.getGammaDimension();
  const auto& matrices = dirac.getMatrices();
  const auto& epsilons = dirac.getEpsilons();

  const cx_mat m_2 = matrices[k] * matrices[k];
  const cx_mat m_3 = matrices[k] * m_2;

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

}// namespace

// Directional matrix derivative of Tr(D^2):
// grad_{M_k} Tr(D^2) = 4 * d_gamma * [ N * M_k + eps_k * Tr(M_k) * I_N ].
arma::cx_mat BarrettGlaserAction::derDirac2(const IDiracOperator& dirac, const int k) {
  const auto dim = dirac.getMatrixDimension();
  const auto gamma_dim = dirac.getGammaDimension();
  const auto& epsilons = dirac.getEpsilons();
  const auto& matrices = dirac.getMatrices();

  arma::cx_mat res(dim, dim, fill::eye);
  res *= epsilons[k] * trace(matrices[k]).real();
  res += dim * matrices[k];
  return 4 * gamma_dim * res;
}

// Directional matrix derivative of Tr(D^4):
// grad_{M_k} Tr(D^4) = B4(k) + B2(k) + B(k).
// When herm is true, projects onto Hermitian tangent space: 2 * (grad + grad†).
arma::cx_mat BarrettGlaserAction::derDirac4(const IDiracOperator& dirac, const int k, const bool herm) {
  const auto dim = dirac.getMatrixDimension();
  const auto num_matrices = dirac.getNumMatrices();
  const auto& epsilons = dirac.getEpsilons();
  const auto& omega_table_4 = dirac.getOmegaTable4();

  cx_mat res(dim, dim, fill::zeros);

  for (int i_1 = 0; i_1 < num_matrices; ++i_1) {
    if (i_1 != k) {
      for (int i_2 = i_1 + 1; i_2 < num_matrices; ++i_2) {
        if (i_2 != k) {
          for (int i_3 = i_2 + 1; i_3 < num_matrices; ++i_3) {
            if (i_3 != k) {
              if (double e = epsilons[k] * epsilons[i_1] * epsilons[i_2] * epsilons[i_3]; e < 0) {
                double cliff_1 = omega_table_4[i_3 + num_matrices * (i_2 + num_matrices * (i_1 + num_matrices * k))].imag();
                double cliff_2 = omega_table_4[i_2 + num_matrices * (i_3 + num_matrices * (i_1 + num_matrices * k))].imag();
                double cliff_3 = omega_table_4[i_3 + num_matrices * (i_1 + num_matrices * (i_2 + num_matrices * k))].imag();

                if (std::abs(cliff_1) > 1e-10) {
                  res += computeB4(dirac, k, i_1, i_2, i_3, cliff_1, true);
                  res += computeB4(dirac, k, i_1, i_3, i_2, cliff_2, true);
                  res += computeB4(dirac, k, i_2, i_1, i_3, cliff_3, true);
                }
              } else {
                double cliff_1 = omega_table_4[i_3 + num_matrices * (i_2 + num_matrices * (i_1 + num_matrices * k))].real();
                double cliff_2 = omega_table_4[i_2 + num_matrices * (i_3 + num_matrices * (i_1 + num_matrices * k))].real();
                double cliff_3 = omega_table_4[i_3 + num_matrices * (i_1 + num_matrices * (i_2 + num_matrices * k))].real();

                if (std::abs(cliff_1) > 1e-10) {
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

  for (int i = 0; i < num_matrices; ++i) {
    if (i != k) {
      res += computeB2(dirac, k, i);
    }
  }

  res += computeB(dirac, k);

  if (herm) {
    // Projects the gradient onto the Hermitian tangent space: grad_H = 2 * (grad + grad†).
    // The factor of 2 accounts for symmetric variations of off-diagonal elements (EP-1).
    return 2 * (res + res.t());
  } else {
    return 4 * res;
  }
}

// Evaluates the full action gradient:
// grad_{M_k} S(D) = g_2 * grad_{M_k} Tr(D^2) + g_4 * grad_{M_k} Tr(D^4).
arma::cx_mat BarrettGlaserAction::derDirac24(
    const IDiracOperator& dirac,
    const int k,
    const bool herm,
    const double g_2,
    const double g_4) {
  return g_2 * derDirac2(dirac, k) + g_4 * derDirac4(dirac, k, herm);
}

arma::cx_mat BarrettGlaserAction::derDirac24(
    const IDiracOperator& dirac,
    const int k,
    const bool herm) const {
  return derDirac24(dirac, k, herm, m_g_2, m_g_4);
}
