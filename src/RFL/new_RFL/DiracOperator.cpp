//
// Created by Paul Druce on 12/11/2022.
//

#include "DiracOperator.hpp"
#include "Clifford.hpp"
#include <cassert>

using namespace std;
using namespace arma;

DiracOperator::DiracOperator(int p, int q, int dim)
    : m_clifford(Clifford(p, q)), m_dim(dim) {
  int n = p + q;

  vector<cx_mat> gamma = m_clifford.getGammaMatrices();
  this->m_gamma_dim = m_clifford.getGammaDimension();

  vector<cx_mat> herm;
  herm.reserve(p);
  for (int i = 0; i < p; ++i) {
    herm.push_back(gamma[i]);
  }

  vector<cx_mat> anti;
  anti.reserve(q);
  for (int i = 0; i < q; ++i) {
    anti.emplace_back(cx_double(0, 1) * gamma[p + i]);
  }

  int count = (int)pow(2, n);
  // The outer for loop will run 2^n times (the number of all possible subsets).
  // Here variable i will act as a binary counter
  for (int i = 0; i < count; i++) {
    vector<int> vec;
    // The inner for loop will run n times, As the maximum number of elements a set can have is n
    // This loop will generate a subset
    for (int j = 0; j < n; j++) {
      // This if condition will check if jth bit in binary representation of i is set or not
      // if the value of (i & (1 << j)) is greater than 0, include arr[j] in the current subset
      // otherwise exclude arr[j]
      if ((i & (1 << j)) > 0) {
        vec.push_back(j);
      }
    }

    // Now calculate and push product if it has odd number of gammas
    int k = (int)vec.size();
    if (k % 2 && k != 1) {
      vector<int>::const_iterator begin(vec.begin());
      //            vector<int>::const_iterator end(vec.end());
      cx_mat mat = gamma.at(*begin);
      //            for (auto iter = vec.begin() + 1; iter != end; ++iter) {
      //                M *= gamma.at((*iter));
      //            }
      bool first = true;
      for (const auto& v : vec) {
        if (first) {// skip first entry
          first = false;
          continue;
        }
        mat *= gamma.at(v);
      }
      if (mat.is_hermitian()) {
        herm.push_back(mat);
      } else {
        anti.emplace_back(cx_double(0, 1) * mat);
      }
    }
  }

  this->m_num_herm = (int)herm.size();
  this->m_num_antiherm = (int)anti.size();
  this->m_num_matrices = m_num_herm + m_num_antiherm;

  this->m_omegas = std::make_unique<std::vector<cx_mat>>(m_num_matrices);
  for (int i = 0; i < m_num_herm; ++i) {
    (*m_omegas)[i] = herm[i];
  }
  for (int i = 0; i < m_num_antiherm; ++i) {
    (*m_omegas)[m_num_herm + i] = anti[i];
  }

  initOmegaTable4();

  // allocate and initialize H and L matrices to identity
  m_matrices = std::make_unique<std::vector<arma::cx_mat>>(m_num_matrices);
  m_momenta = std::make_unique<std::vector<arma::cx_mat>>(m_num_matrices);
  m_epsilons = std::make_unique<std::vector<int>>(m_num_matrices);
  for (int i = 0; i < m_num_matrices; i++) {
    if (i < m_num_herm) {
      (*m_epsilons)[i] = 1;
    } else {
      (*m_epsilons)[i] = -1;
    }

    (*m_matrices)[i].eye(m_dim, m_dim);
    (*m_momenta)[i].eye(m_dim, m_dim);
  }
}

DiracOperator::DiracOperator(const DiracOperator& original)
    : m_clifford(original.m_clifford) {
  // Make a copy of all state of original DiracOperator.
  this->m_dim = original.m_dim;
  this->m_num_matrices = original.m_num_herm;
  this->m_num_herm = original.m_num_herm;
  this->m_num_antiherm = original.m_num_antiherm;
  this->m_gamma_dim = original.m_gamma_dim;
  this->m_matrices = std::make_unique<std::vector<arma::cx_mat>>(*original.m_matrices);
  this->m_momenta = std::make_unique<std::vector<arma::cx_mat>>(*original.m_momenta);
  this->m_omegas = std::make_unique<std::vector<cx_mat>>(*original.m_omegas);
  this->m_epsilons = std::make_unique<std::vector<int>>(*original.m_epsilons);
  this->m_omega_table_4 = std::make_unique<std::vector<cx_double>>(*original.m_omega_table_4);
}

// TODO: figure out what this function is for.
/**
 * @brief standalone function to do some form of conversion.
 *
 * @param dec
 * @param base
 * @param max
 * @return vector<int>
 */
static vector<int> baseConversion(int dec, const int& base, const int& max) {
  vector<int> rem;

  while (dec) {
    rem.push_back(dec % base);
    dec /= base;
  }

  for (int i = (int)rem.size(); i < max; i++)
    rem.push_back(0);

  reverse(rem.begin(), rem.end());

  return rem;
}

cx_mat DiracOperator::getDiracMatrix() const {
  // initialize dirac op to zero
  const int dim_dirac = m_dim * m_dim * m_gamma_dim;
  cx_mat dirac(dim_dirac, dim_dirac, fill::zeros);

  const cx_mat id(m_dim, m_dim, fill::eye);
  for (int i = 0; i < m_num_matrices; ++i) {
    cx_mat bracket = kron((*m_matrices)[i], id) + (*m_epsilons)[i] * kron(id, (*m_matrices)[i].st());
    dirac += kron((*m_omegas)[i], bracket);
  }

  return dirac;
}

vec DiracOperator::getEigenvalues() const {
  const auto diracMatrix = getDiracMatrix();
  assert(diracMatrix.is_hermitian(1e-16));
  auto eigen_vals = arma::eig_sym(diracMatrix);
  return eigen_vals;
}

vector<cx_mat> DiracOperator::getHermitianMatrices() const {
  std::vector<arma::cx_mat> herm_matrices;
  for (std::size_t i = 0; i < m_matrices->size(); i++) {
    if (m_epsilons->at(i) == 1) {
      herm_matrices.push_back(m_matrices->at(i));
    }
  }
  return herm_matrices;
}

vector<cx_mat> DiracOperator::getAntiHermitianMatrices() const {
  vector<arma::cx_mat> anti_herm_matrices;
  for (size_t i = 0; i < m_matrices->size(); i++) {
    if (m_epsilons->at(i) == -1) {
      anti_herm_matrices.emplace_back(m_matrices->at(i) * cx_double(0, 1));
    }
  }
  return anti_herm_matrices;
}

void DiracOperator::printOmegaTable4() const {
  const int n = (int)pow(m_num_matrices, 4);

  for (int i = 0; i < n; ++i) {
    cx_double z = (*m_omega_table_4)[i];
    if (z != cx_double(0., 0.)) {
      int e = 1;
      vector<int> prod = baseConversion(i, m_num_matrices, 4);
      //            vector<int>::const_iterator end(prod.end());
      //            for (vector<int>::const_iterator iter = prod.begin(); iter != end; ++iter) {
      //                cout << (*iter) << " ";
      //                e *= eps[(*iter)];
      //            }
      for (const auto& p : prod) {
        cout << p << " ";
        e *= (*m_epsilons)[p];
      }
      cout << " " << (*m_omega_table_4)[i] << e << endl;
    }
  }
}

void DiracOperator::initOmegaTable4() {
  this->m_omega_table_4 = std::make_unique<std::vector<cx_double>>(m_num_matrices * m_num_matrices * m_num_matrices * m_num_matrices);

  for (int i = 0; i < m_num_matrices; ++i) {
    for (int j = 0; j < m_num_matrices; ++j) {
      for (int k = 0; k < m_num_matrices; ++k) {
        for (int l = 0; l < m_num_matrices; ++l)
          (*m_omega_table_4)[l + m_num_matrices * (k + m_num_matrices * (j + m_num_matrices * i))] = trace(
              (*m_omegas)[i] * (*m_omegas)[j] * (*m_omegas)[k] * (*m_omegas)[l]);
      }
    }
  }
}

cx_mat DiracOperator::derDirac24(const int& k, const bool& herm, const double g_2) const {
  return g_2 * derDirac2(k) + derDirac4(k, herm);
}

cx_mat DiracOperator::derDirac2(const int& k) const {
  cx_mat res(m_dim, m_dim, fill::eye);

  res *= (*m_epsilons)[k] * trace((*m_matrices)[k]).real();
  res += m_dim * (*m_matrices)[k];

  return 4 * m_gamma_dim * res;
}

cx_mat DiracOperator::derDirac4(const int& k, const bool& herm) const {
  cx_mat res(m_dim, m_dim, fill::zeros);
  auto& epsilons = this->getEpsilons();
  auto& omega_table_4 = this->getOmegaTable4();

  // four distinct indices
  for (int i_1 = 0; i_1 < m_num_matrices; ++i_1) {
    if (i_1 != k) {
      for (int i_2 = i_1 + 1; i_2 < m_num_matrices; ++i_2) {
        if (i_2 != k) {
          for (int i_3 = i_2 + 1; i_3 < m_num_matrices; ++i_3) {
            if (i_3 != k) {
              // epsilon factor
              double e = epsilons[k] * epsilons[i_1] * epsilons[i_2] * epsilons[i_3];

              if (e < 0) {
                // clifford product
                double cliff_1 = omega_table_4[i_3 + m_num_matrices * (i_2 + m_num_matrices * (i_1 + m_num_matrices * k))].imag();
                double cliff_2 = omega_table_4[i_2 + m_num_matrices * (i_3 + m_num_matrices * (i_1 + m_num_matrices * k))].imag();
                double cliff_3 = omega_table_4[i_3 + m_num_matrices * (i_1 + m_num_matrices * (i_2 + m_num_matrices * k))].imag();

                if (fabs(cliff_1) > 1e-10) {
                  res += computeB4(k, i_1, i_2, i_3, cliff_1, true);
                  res += computeB4(k, i_1, i_3, i_2, cliff_2, true);
                  res += computeB4(k, i_2, i_1, i_3, cliff_3, true);
                }
              } else {
                // clifford product
                double cliff_1 = omega_table_4[i_3 + m_num_matrices * (i_2 + m_num_matrices * (i_1 + m_num_matrices * k))].real();
                double cliff_2 = omega_table_4[i_2 + m_num_matrices * (i_3 + m_num_matrices * (i_1 + m_num_matrices * k))].real();
                double cliff_3 = omega_table_4[i_3 + m_num_matrices * (i_1 + m_num_matrices * (i_2 + m_num_matrices * k))].real();

                if (fabs(cliff_1) > 1e-10) {
                  res += computeB4(k, i_1, i_2, i_3, cliff_1, false);
                  res += computeB4(k, i_1, i_3, i_2, cliff_2, false);
                  res += computeB4(k, i_2, i_1, i_3, cliff_3, false);
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
  for (int i = 0; i < m_num_matrices; ++i) {
    if (i != k) {
      res += computeB2(k, i);
    }
  }

  // all indices equal
  res += computeB(k);

  cx_mat final_result;
  if (herm) {
    final_result = 2 * (res + res.t());
  } else {
    final_result = 4 * res;
  }

  return final_result;
}

cx_mat DiracOperator::computeB4(const int& k,
                                const int& i_2,
                                const int& i_3,
                                const int& i_4,
                                const double& cliff,
                                const bool& neg) const {
  auto& matrices = this->getMatrices();
  auto& epsilons = this->getEpsilons();

  // base matrix products
  cx_mat m_2_m_3 = matrices[i_2] * matrices[i_3];
  cx_mat m_2_m_4 = matrices[i_2] * matrices[i_4];
  cx_mat m_3_m_4 = matrices[i_3] * matrices[i_4];
  cx_mat m_2_m_3_m_4 = m_2_m_3 * matrices[i_4];

  // return value
  cx_mat res(m_dim, m_dim, fill::eye);

  if (neg) {

    // traces
    double tr_234 = trace(m_2_m_3_m_4).imag();
    double tr_2 = trace(matrices[i_2]).real();
    double tr_3 = trace(matrices[i_3]).real();
    double tr_4 = trace(matrices[i_4]).real();

    // compute sum
    res *= -2 * epsilons[k] * tr_234;
    res += cx_double(0., m_dim) * (m_2_m_3_m_4 - m_2_m_3_m_4.t());
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
    res += m_dim * (m_2_m_3_m_4 + m_2_m_3_m_4.t());
    res += epsilons[i_2] * tr_2 * (m_3_m_4 + m_3_m_4.t());
    res += epsilons[i_3] * tr_3 * (m_2_m_4 + m_2_m_4.t());
    res += epsilons[i_4] * tr_4 * (m_2_m_3 + m_2_m_3.t());
    res += 2 * epsilons[k] * epsilons[i_2] * tr_34 * matrices[i_2];
    res += 2 * epsilons[k] * epsilons[i_3] * tr_24 * matrices[i_3];
    res += 2 * epsilons[k] * epsilons[i_4] * tr_23 * matrices[i_4];
  }

  return cliff * res;
}

cx_mat DiracOperator::computeB2(const int& k, const int& i) const {
  auto& omega_table_4 = this->getOmegaTable4();
  auto& matrices = this->getMatrices();
  auto& epsilons = this->getEpsilons();

  // clifford product
  double cliff = omega_table_4[i + m_num_matrices * (k + m_num_matrices * (i + m_num_matrices * k))].real();

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
  cx_mat res(m_dim, m_dim, fill::eye);

  if (cliff < 0) {
    // compute sum
    res *= epsilons[k] * triki;
    res += m_dim * (mi_mi_mk + mi_mi_mk.t() - mi_mk_mi);
    res += epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 2 * epsilons[k] * epsilons[i] * trik * matrices[i];
    res += epsilons[k] * trk * mi_mi;
    res += trii * matrices[k];
  } else {
    // compute sum
    res *= 3 * epsilons[k] * triki;
    res += m_dim * (mi_mi_mk + mi_mi_mk.t() + mi_mk_mi);
    res += 3 * epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 6 * epsilons[k] * epsilons[i] * trik * matrices[i];
    res += 3 * epsilons[k] * trk * mi_mi;
    res += 3 * trii * matrices[k];
  }

  return 2 * m_gamma_dim * res;
}

cx_mat DiracOperator::computeB(const int& k) const {
  const auto& matrices = this->getMatrices();
  const auto& epsilons = this->getEpsilons();

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

// TODO: Refactor this method, to reuse the randomisation process. Can just pass in the matrix or array of matrices to assign
void DiracOperator::randomiseMatrices(const IRng& rng) const {
  auto& matrices = this->getMatrices();
  auto& momenta = this->getMomenta();

  for (int i = 0; i < m_num_matrices; ++i) {
    // loop on indices
    for (int j = 0; j < m_dim; ++j) {
      const double x = rng.getGaussian(1.0);
      matrices[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < m_dim; ++k) {
        const double a = rng.getGaussian(1.0);
        const double b = rng.getGaussian(1.0);
        matrices[i](j, k) = cx_double(a, b) / sqrt(2.);
        matrices[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }

  for (int i = 0; i < m_num_matrices; ++i) {
    // loop on indices
    for (int j = 0; j < m_dim; ++j) {
      const double x = rng.getGaussian(1.0);
      momenta[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < m_dim; ++k) {
        const double a = rng.getGaussian(1.0);
        const double b = rng.getGaussian(1.0);
        momenta[i](j, k) = cx_double(a, b) / sqrt(2.);
        momenta[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }
}

double DiracOperator::traceOfDiracSquared() const {
  const auto& matrices = this->getMatrices();
  const auto& epsilons = this->getEpsilons();

  double res = 0.;
  for (int i = 0; i < m_num_matrices; ++i) {
    const double trace_squared = trace(matrices[i] * matrices[i]).real();
    const double trace_mat = trace(matrices[i]).real();

    res += (m_dim * trace_squared + epsilons[i] * trace_mat * trace_mat);
  }

  return 2. * m_gamma_dim * res;
}

double DiracOperator::traceOfDirac4() const {
  double res = 0.;

  const auto num_matrices = this->getNumMatrices();

  // four distinct indices
  for (int i = 0; i < num_matrices; ++i) {
    for (int j = i + 1; j < num_matrices; ++j) {
      for (int k = j + 1; k < num_matrices; ++k) {
        for (int l = k + 1; l < num_matrices; ++l)
          res += 8 * (computeA4(i, j, k, l) + computeA4(i, j, l, k) + computeA4(i, k, j, l));
      }
    }
  }

  // two distinct pairs of equal indices
  for (int i = 0; i < num_matrices; ++i) {
    for (int j = i + 1; j < num_matrices; ++j)
      res += 2 * computeA2(i, j);
  }

  // all indices equal
  for (int i = 0; i < num_matrices; ++i)
    res += computeA(i);

  return res;
}

double DiracOperator::computeA4(const int& i_1, const int& i_2, const int& i_3, const int& i_4) const {
  auto& epsilons = this->getEpsilons();
  auto& omega_table_4 = this->getOmegaTable4();
  auto& matrices = this->getMatrices();

  const int e = epsilons[i_1] * epsilons[i_2] * epsilons[i_3] * epsilons[i_4];

  // if e=-1, then [1+*e] becomes 2i*imag
  // and the clifford part is guaranteed to
  // be pure imaginary
  if (e < 0) {
    // clifford product
    const double cliff = omega_table_4[i_4 + m_num_matrices * (i_3 + m_num_matrices * (i_2 + m_num_matrices * i_1))].imag();

    if (fabs(cliff) > 1e-10) {
      // base matrix products
      const cx_mat m_1_m_2 = matrices[i_1] * matrices[i_2];
      const cx_mat m_1_m_3 = matrices[i_1] * matrices[i_3];
      const cx_mat m_1_m_4 = matrices[i_1] * matrices[i_4];
      const cx_mat m_2_m_3 = matrices[i_2] * matrices[i_3];
      const cx_mat m_2_m_4 = matrices[i_2] * matrices[i_4];
      const cx_mat m_3_m_4 = matrices[i_3] * matrices[i_4];

      // traces
      double tr_1234 = trace(m_1_m_2 * m_3_m_4).imag();
      double tr_234 = trace(m_2_m_3 * matrices[i_4]).imag();
      double tr_134 = trace(m_1_m_3 * matrices[i_4]).imag();
      double tr_124 = trace(m_1_m_2 * matrices[i_4]).imag();
      double tr_123 = trace(m_1_m_2 * matrices[i_3]).imag();
      double tr_1 = trace(matrices[i_1]).real();
      double tr_2 = trace(matrices[i_2]).real();
      double tr_3 = trace(matrices[i_3]).real();
      double tr_4 = trace(matrices[i_4]).real();

      // compute sum
      double res = m_dim * tr_1234;
      res += epsilons[i_1] * tr_1 * tr_234;
      res += epsilons[i_2] * tr_2 * tr_134;
      res += epsilons[i_3] * tr_3 * tr_124;
      res += epsilons[i_4] * tr_4 * tr_123;

      return -2 * cliff * res;
      // NOTE: this minus here comes from the 'i' in cliff multiplied by
      // the 'i' coming from 2i*imag.
    } else {
      return 0.;
    }
  } else {
    // clifford product
    const double cliff = omega_table_4[i_4 + m_num_matrices * (i_3 + m_num_matrices * (i_2 + m_num_matrices * i_1))].real();

    if (fabs(cliff) > 1e-10) {
      // base matrix products
      cx_mat m_1_m_2 = matrices[i_1] * matrices[i_2];
      cx_mat m_1_m_3 = matrices[i_1] * matrices[i_3];
      cx_mat m_1_m_4 = matrices[i_1] * matrices[i_4];
      cx_mat m_2_m_3 = matrices[i_2] * matrices[i_3];
      cx_mat m_2_m_4 = matrices[i_2] * matrices[i_4];
      cx_mat m_3_m_4 = matrices[i_3] * matrices[i_4];

      // traces
      double tr_1234 = trace(m_1_m_2 * m_3_m_4).real();
      double tr_234 = trace(m_2_m_3 * matrices[i_4]).real();
      double tr_134 = trace(m_1_m_3 * matrices[i_4]).real();
      double tr_124 = trace(m_1_m_2 * matrices[i_4]).real();
      double tr_123 = trace(m_1_m_2 * matrices[i_3]).real();
      double tr_12 = trace(m_1_m_2).real();
      double tr_34 = trace(m_3_m_4).real();
      double tr_13 = trace(m_1_m_3).real();
      double tr_24 = trace(m_2_m_4).real();
      double tr_14 = trace(m_1_m_4).real();
      double tr_23 = trace(m_2_m_3).real();
      double tr_1 = trace(matrices[i_1]).real();
      double tr_2 = trace(matrices[i_2]).real();
      double tr_3 = trace(matrices[i_3]).real();
      double tr_4 = trace(matrices[i_4]).real();

      double res = m_dim * tr_1234;
      res += epsilons[i_1] * tr_1 * tr_234;
      res += epsilons[i_2] * tr_2 * tr_134;
      res += epsilons[i_3] * tr_3 * tr_124;
      res += epsilons[i_4] * tr_4 * tr_123;
      res += epsilons[i_1] * epsilons[i_2] * tr_12 * tr_34;
      res += epsilons[i_1] * epsilons[i_3] * tr_13 * tr_24;
      res += epsilons[i_1] * epsilons[i_4] * tr_14 * tr_23;

      //	  cliff = omega_table_4[i4 + D.nHL * (i3 + D.nHL * (i2 + D.nHL * i1))].real();

      return 2 * cliff * res;
    } else {
      return 0.;
    }
  }
}

double DiracOperator::computeA2(const int& i_1, const int& i_2) const {
  const auto& omega_table_4 = this->getOmegaTable4();
  const auto& matrices = this->getMatrices();
  const auto& epsilons = this->getEpsilons();

  // clifford product
  const double cliff = omega_table_4[i_2 + m_num_matrices * (i_1 + m_num_matrices * (i_2 + m_num_matrices * i_1))].real();

  // base matrix products
  const cx_mat m_1_m_1 = matrices[i_1] * matrices[i_1];
  const cx_mat m_2_m_2 = matrices[i_2] * matrices[i_2];
  const cx_mat m_1_m_2 = matrices[i_1] * matrices[i_2];

  // traces
  const double tr_1122 = trace(m_1_m_1 * m_2_m_2).real();
  const double tr_1212 = trace(m_1_m_2 * m_1_m_2).real();
  const double tr_122 = trace(m_1_m_2 * matrices[i_2]).real();
  const double tr_112 = trace(m_1_m_1 * matrices[i_2]).real();
  const double tr_11 = trace(m_1_m_1).real();
  const double tr_22 = trace(m_2_m_2).real();
  const double tr_12 = trace(m_1_m_2).real();
  const double tr_1 = trace(matrices[i_1]).real();
  const double tr_2 = trace(matrices[i_2]).real();

  if (cliff < 0) {
    // compute sum
    double res = m_dim * (2 * tr_1122 - tr_1212);
    res += 2 * epsilons[i_1] * tr_1 * tr_122;
    res += 2 * epsilons[i_2] * tr_2 * tr_112;
    res += tr_11 * tr_22;
    res += 2 * epsilons[i_1] * epsilons[i_2] * tr_12 * tr_12;

    return 2 * m_gamma_dim * res;
  } else {
    // compute sum
    double res = m_dim * (2 * tr_1122 + tr_1212);
    res += 6 * epsilons[i_1] * tr_1 * tr_122;
    res += 6 * epsilons[i_2] * tr_2 * tr_112;
    res += 3 * tr_11 * tr_22;
    res += 6 * epsilons[i_1] * epsilons[i_2] * tr_12 * tr_12;

    return 2 * m_gamma_dim * res;
  }
}

double DiracOperator::computeA(const int& i) const {
  const auto& matrices = this->getMatrices();
  const auto& epsilons = this->getEpsilons();

  // base matrix products
  const cx_mat m_2 = matrices[i] * matrices[i];
  const cx_mat m_3 = m_2 * matrices[i];

  // traces
  const double tr_1 = trace(matrices[i]).real();
  const double tr_2 = trace(m_2).real();
  const double tr_3 = trace(m_3).real();
  const double tr_4 = trace(m_3 * matrices[i]).real();

  double res = m_dim * tr_4;
  res += 4 * epsilons[i] * tr_1 * tr_3;
  res += 3 * tr_2 * tr_2;

  return 2 * m_gamma_dim * res;
}
