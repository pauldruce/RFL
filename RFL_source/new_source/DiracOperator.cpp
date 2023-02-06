//
// Created by Paul Druce on 12/11/2022.
//

#include "DiracOperator.hpp"
#include "Clifford.hpp"

using namespace std;
using namespace arma;

DiracOperator::DiracOperator(int p, int q, int dim)
    : m_dim(dim) {
  int n = p + q;

  // create a type (p, q) clifford module
  Clifford clifford(p, q);
  vector<cx_mat> gamma = clifford.getGammaMatrices();
  this->m_gamma_dim = clifford.getGammaDimension();

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
        if (first) {// skipp first entry
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

  this->m_omegas = new cx_mat[m_num_matrices];
  for (int i = 0; i < m_num_herm; ++i) {
    m_omegas[i] = herm[i];
  }
  for (int i = 0; i < m_num_antiherm; ++i) {
    m_omegas[m_num_herm + i] = anti[i];
  }

  initOmegaTable4();

  // allocate and initialize H and L matrices to identity
  m_matrices = new arma::cx_mat[m_num_matrices];
  m_momenta = new arma::cx_mat[m_num_matrices];
  m_epsilons = new int[m_num_matrices];
  for (int i = 0; i < m_num_matrices; i++) {
    if (i < m_num_herm) {
      m_epsilons[i] = 1;
    } else {
      m_epsilons[i] = -1;
    }

    m_matrices[i].eye(m_dim, m_dim);
    m_momenta[i].eye(m_dim, m_dim);
  }
}

DiracOperator::~DiracOperator() {
  delete[] m_matrices;
  delete[] m_momenta;
  delete[] m_epsilons;
  delete[] m_omegas;
  delete[] m_omega_table_4;
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
vector<int> baseConversion(int dec, const int& base, const int& max) {
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
  int dim_dirac = m_dim * m_dim * m_gamma_dim;
  cx_mat dirac(dim_dirac, dim_dirac, fill::zeros);

  cx_mat id(m_dim, m_dim, fill::eye);
  for (int i = 0; i < m_num_matrices; ++i) {
    cx_mat bracket = kron(m_matrices[i], id) + m_epsilons[i] * kron(id, m_matrices[i].st());
    dirac += kron(m_omegas[i], bracket);
  }

  return dirac;
}

void DiracOperator::printOmegaTable4() const {
  const int n = (int)pow(m_num_matrices, 4);

  for (int i = 0; i < n; ++i) {
    cx_double z = m_omega_table_4[i];
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
        e *= m_epsilons[p];
      }
      cout << " " << m_omega_table_4[i] << e << endl;
    }
  }
}

void DiracOperator::initOmegaTable4() {
  this->m_omega_table_4 = new cx_double[m_num_matrices * m_num_matrices * m_num_matrices * m_num_matrices];

  for (int i = 0; i < m_num_matrices; ++i) {
    for (int j = 0; j < m_num_matrices; ++j) {
      for (int k = 0; k < m_num_matrices; ++k) {
        for (int l = 0; l < m_num_matrices; ++l)
          m_omega_table_4[l + m_num_matrices * (k + m_num_matrices * (j + m_num_matrices * i))] = trace(
              m_omegas[i] * m_omegas[j] * m_omegas[k] * m_omegas[l]);
      }
    }
  }
}

cx_mat DiracOperator::derDirac24(const int& k, const bool& herm, const double g_2) const {
  return g_2 * derDirac2(k) + derDirac4(k, herm);
}

cx_mat DiracOperator::derDirac2(const int& k) const {
  cx_mat res(m_dim, m_dim, fill::eye);

  res *= m_epsilons[k] * trace(m_matrices[k]).real();
  res += m_dim * m_matrices[k];

  return 4 * m_gamma_dim * res;
}

cx_mat DiracOperator::derDirac4(const int& k, const bool& herm) const {
  cx_mat res(m_dim, m_dim, fill::zeros);

  // four distinct indices
  for (int i_1 = 0; i_1 < m_num_matrices; ++i_1) {
    if (i_1 != k) {
      for (int i_2 = i_1 + 1; i_2 < m_num_matrices; ++i_2) {
        if (i_2 != k) {
          for (int i_3 = i_2 + 1; i_3 < m_num_matrices; ++i_3) {
            if (i_3 != k) {
              // epsilon factor
              double e = m_epsilons[k] * m_epsilons[i_1] * m_epsilons[i_2] * m_epsilons[i_3];

              if (e < 0) {
                // clifford product
                double cliff_1 = m_omega_table_4[i_3 + m_num_matrices * (i_2 + m_num_matrices * (i_1 + m_num_matrices * k))].imag();
                double cliff_2 = m_omega_table_4[i_2 + m_num_matrices * (i_3 + m_num_matrices * (i_1 + m_num_matrices * k))].imag();
                double cliff_3 = m_omega_table_4[i_3 + m_num_matrices * (i_1 + m_num_matrices * (i_2 + m_num_matrices * k))].imag();

                if (fabs(cliff_1) > 1e-10) {
                  res += computeB4(k, i_1, i_2, i_3, cliff_1, true);
                  res += computeB4(k, i_1, i_3, i_2, cliff_2, true);
                  res += computeB4(k, i_2, i_1, i_3, cliff_3, true);
                }
              } else {
                // clifford product
                double cliff_1 = m_omega_table_4[i_3 + m_num_matrices * (i_2 + m_num_matrices * (i_1 + m_num_matrices * k))].real();
                double cliff_2 = m_omega_table_4[i_2 + m_num_matrices * (i_3 + m_num_matrices * (i_1 + m_num_matrices * k))].real();
                double cliff_3 = m_omega_table_4[i_3 + m_num_matrices * (i_1 + m_num_matrices * (i_2 + m_num_matrices * k))].real();

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

  if (herm) {
    return 2 * (res + res.t());
  } else {
    return 4 * res;
  }
}

cx_mat DiracOperator::computeB4(const int& k,
                                const int& i_2,
                                const int& i_3,
                                const int& i_4,
                                const double& cliff,
                                const bool& neg) const {
  // base matrix products
  cx_mat m_2_m_3 = m_matrices[i_2] * m_matrices[i_3];
  cx_mat m_2_m_4 = m_matrices[i_2] * m_matrices[i_4];
  cx_mat m_3_m_4 = m_matrices[i_3] * m_matrices[i_4];
  cx_mat m_2_m_3_m_4 = m_2_m_3 * m_matrices[i_4];

  // return value
  cx_mat res(m_dim, m_dim, fill::eye);

  if (neg) {

    // traces
    double tr_234 = trace(m_2_m_3_m_4).imag();
    double tr_2 = trace(m_matrices[i_2]).real();
    double tr_3 = trace(m_matrices[i_3]).real();
    double tr_4 = trace(m_matrices[i_4]).real();

    // compute sum
    res *= -2 * m_epsilons[k] * tr_234;
    res += cx_double(0., m_dim) * (m_2_m_3_m_4 - m_2_m_3_m_4.t());
    res += cx_double(0., m_epsilons[i_2] * tr_2) * (m_3_m_4 - m_3_m_4.t());
    res += cx_double(0., m_epsilons[i_3] * tr_3) * (m_2_m_4 - m_2_m_4.t());
    res += cx_double(0., m_epsilons[i_4] * tr_4) * (m_2_m_3 - m_2_m_3.t());
  } else {
    // traces
    double tr_234 = trace(m_2_m_3_m_4).real();
    double tr_23 = trace(m_2_m_3).real();
    double tr_24 = trace(m_2_m_4).real();
    double tr_34 = trace(m_3_m_4).real();
    double tr_2 = trace(m_matrices[i_2]).real();
    double tr_3 = trace(m_matrices[i_3]).real();
    double tr_4 = trace(m_matrices[i_4]).real();

    // compute sum
    res *= 2 * m_epsilons[k] * tr_234;
    res += m_dim * (m_2_m_3_m_4 + m_2_m_3_m_4.t());
    res += m_epsilons[i_2] * tr_2 * (m_3_m_4 + m_3_m_4.t());
    res += m_epsilons[i_3] * tr_3 * (m_2_m_4 + m_2_m_4.t());
    res += m_epsilons[i_4] * tr_4 * (m_2_m_3 + m_2_m_3.t());
    res += 2 * m_epsilons[k] * m_epsilons[i_2] * tr_34 * m_matrices[i_2];
    res += 2 * m_epsilons[k] * m_epsilons[i_3] * tr_24 * m_matrices[i_3];
    res += 2 * m_epsilons[k] * m_epsilons[i_4] * tr_23 * m_matrices[i_4];
  }

  return cliff * res;
}

cx_mat DiracOperator::computeB2(const int& k, const int& i) const {
  // clifford product
  double cliff = m_omega_table_4[i + m_num_matrices * (k + m_num_matrices * (i + m_num_matrices * k))].real();

  // base matrix products
  cx_mat mi_mk = m_matrices[i] * m_matrices[k];
  cx_mat mi_mi = m_matrices[i] * m_matrices[i];
  cx_mat mi_mi_mk = m_matrices[i] * mi_mk;
  cx_mat mi_mk_mi = mi_mk * m_matrices[i];

  // traces
  double triki = trace(mi_mk_mi).real();
  double trik = trace(mi_mk).real();
  double trii = trace(mi_mi).real();
  double tri = trace(m_matrices[i]).real();
  double trk = trace(m_matrices[k]).real();

  // return value
  cx_mat res(m_dim, m_dim, fill::eye);

  if (cliff < 0) {
    // compute sum
    res *= m_epsilons[k] * triki;
    res += m_dim * (mi_mi_mk + mi_mi_mk.t() - mi_mk_mi);
    res += m_epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 2 * m_epsilons[k] * m_epsilons[i] * trik * m_matrices[i];
    res += m_epsilons[k] * trk * mi_mi;
    res += trii * m_matrices[k];
  } else {
    // compute sum
    res *= 3 * m_epsilons[k] * triki;
    res += m_dim * (mi_mi_mk + mi_mi_mk.t() + mi_mk_mi);
    res += 3 * m_epsilons[i] * tri * (mi_mk + mi_mk.t());
    res += 6 * m_epsilons[k] * m_epsilons[i] * trik * m_matrices[i];
    res += 3 * m_epsilons[k] * trk * mi_mi;
    res += 3 * trii * m_matrices[k];
  }

  return 2 * m_gamma_dim * res;
}

cx_mat DiracOperator::computeB(const int& k) const {
  // base matrix products
  cx_mat m_2 = m_matrices[k] * m_matrices[k];
  cx_mat m_3 = m_matrices[k] * m_2;

  // traces
  double tr_3 = trace(m_3).real();
  double tr_2 = trace(m_2).real();
  double tr_1 = trace(m_matrices[k]).real();

  cx_mat res(m_dim, m_dim, fill::eye);
  res *= m_epsilons[k] * tr_3;
  res += m_dim * m_3;
  res += 3 * tr_2 * m_matrices[k];
  res += 3 * m_epsilons[k] * tr_1 * m_2;

  return 2 * m_gamma_dim * res;
}

void DiracOperator::randomiseMatrices(gsl_rng* engine) {
  for (int i = 0; i < m_num_matrices; ++i) {
    // loop on indices
    for (int j = 0; j < m_dim; ++j) {
      double x;
      x = gsl_ran_gaussian(engine, 1.);
      m_matrices[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < m_dim; ++k) {
        double a, b;
        a = gsl_ran_gaussian(engine, 1.);
        b = gsl_ran_gaussian(engine, 1.);
        m_matrices[i](j, k) = cx_double(a, b) / sqrt(2.);
        m_matrices[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }

  for (int i = 0; i < m_num_matrices; ++i) {
    // loop on indices
    for (int j = 0; j < m_dim; ++j) {
      double x;
      x = gsl_ran_gaussian(engine, 1.);
      m_momenta[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < m_dim; ++k) {
        double a, b;
        a = gsl_ran_gaussian(engine, 1.);
        b = gsl_ran_gaussian(engine, 1.);
        m_momenta[i](j, k) = cx_double(a, b) / sqrt(2.);
        m_momenta[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }
}
