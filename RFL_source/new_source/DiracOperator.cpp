//
// Created by Paul Druce on 12/11/2022.
//

#include "DiracOperator.hpp"
#include "Clifford.hpp"

using namespace std;
using namespace arma;

DiracOperator::DiracOperator(int p, int q, int dim_)
    : dim(dim_) {
  int n = p + q;

  // create a type (p, q) clifford module
  Clifford C(p, q);
  vector<cx_mat> gamma = C.getGammaMatrices();
  this->dim_omega = C.getGammaDimension();

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
      cx_mat M = gamma.at(*begin);
      //            for (auto iter = vec.begin() + 1; iter != end; ++iter) {
      //                M *= gamma.at((*iter));
      //            }
      bool first = true;
      for (const auto& v : vec) {
        if (first) {// skipp first entry
          first = false;
          continue;
        }
        M *= gamma.at(v);
      }
      if (M.is_hermitian()) {
        herm.push_back(M);
      } else {
        anti.emplace_back(cx_double(0, 1) * M);
      }
    }
  }

  this->nH = (int)herm.size();
  this->nL = (int)anti.size();
  this->nHL = nH + nL;

  this->omega = new cx_mat[nHL];
  for (int i = 0; i < nH; ++i) {
    omega[i] = herm[i];
  }
  for (int i = 0; i < nL; ++i) {
    omega[nH + i] = anti[i];
  }

  init_omega_table_4();

  // allocate and initialize H and L matrices to identity
  mat = new arma::cx_mat[nHL];
  mom = new arma::cx_mat[nHL];
  eps = new int[nHL];
  for (int i = 0; i < nHL; i++) {
    if (i < nH) {
      eps[i] = 1;
    } else {
      eps[i] = -1;
    }

    mat[i].eye(dim, dim);
    mom[i].eye(dim, dim);
  }
}

DiracOperator::~DiracOperator() {
  delete[] mat;
  delete[] mom;
  delete[] eps;
  delete[] omega;
  delete[] omega_table_4;
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
vector<int> base_conversion(int dec, const int& base, const int& max) {
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

cx_mat DiracOperator::build_dirac() const {
  // initialize dirac op to zero
  int dim_dirac = dim * dim * dim_omega;
  cx_mat dirac(dim_dirac, dim_dirac, fill::zeros);

  cx_mat id(dim, dim, fill::eye);
  for (int i = 0; i < nHL; ++i) {
    cx_mat bracket = kron(mat[i], id) + eps[i] * kron(id, mat[i].st());
    dirac += kron(omega[i], bracket);
  }

  return dirac;
}

void DiracOperator::print_omega_table_4() const {
  const int n = (int)pow(nHL, 4);

  for (int i = 0; i < n; ++i) {
    cx_double z = omega_table_4[i];
    if (z != cx_double(0., 0.)) {
      int e = 1;
      vector<int> prod = base_conversion(i, nHL, 4);
      //            vector<int>::const_iterator end(prod.end());
      //            for (vector<int>::const_iterator iter = prod.begin(); iter != end; ++iter) {
      //                cout << (*iter) << " ";
      //                e *= eps[(*iter)];
      //            }
      for (const auto& p : prod) {
        cout << p << " ";
        e *= eps[p];
      }
      cout << " " << omega_table_4[i] << e << endl;
    }
  }
}

void DiracOperator::init_omega_table_4() {
  this->omega_table_4 = new cx_double[nHL * nHL * nHL * nHL];

  for (int i = 0; i < nHL; ++i) {
    for (int j = 0; j < nHL; ++j) {
      for (int k = 0; k < nHL; ++k) {
        for (int l = 0; l < nHL; ++l)
          omega_table_4[l + nHL * (k + nHL * (j + nHL * i))] = trace(
              omega[i] * omega[j] * omega[k] * omega[l]);
      }
    }
  }
}

cx_mat DiracOperator::der_dirac24(const int& k, const bool& herm, const double g2) const {
  return g2 * der_dirac2(k) + der_dirac4(k, herm);
}

cx_mat DiracOperator::der_dirac2(const int& k) const {
  cx_mat res(dim, dim, fill::eye);

  res *= eps[k] * trace(mat[k]).real();
  res += dim * mat[k];

  return 4 * dim_omega * res;
}

cx_mat DiracOperator::der_dirac4(const int& k, const bool& herm) const {
  cx_mat res(dim, dim, fill::zeros);

  // four distinct indices
  for (int i1 = 0; i1 < nHL; ++i1) {
    if (i1 != k) {
      for (int i2 = i1 + 1; i2 < nHL; ++i2) {
        if (i2 != k) {
          for (int i3 = i2 + 1; i3 < nHL; ++i3) {
            if (i3 != k) {
              // epsilon factor
              double e = eps[k] * eps[i1] * eps[i2] * eps[i3];

              if (e < 0) {
                // clifford product
                double cliff1 = omega_table_4[i3 + nHL * (i2 + nHL * (i1 + nHL * k))].imag();
                double cliff2 = omega_table_4[i2 + nHL * (i3 + nHL * (i1 + nHL * k))].imag();
                double cliff3 = omega_table_4[i3 + nHL * (i1 + nHL * (i2 + nHL * k))].imag();

                if (fabs(cliff1) > 1e-10) {
                  res += compute_B4(k, i1, i2, i3, cliff1, true);
                  res += compute_B4(k, i1, i3, i2, cliff2, true);
                  res += compute_B4(k, i2, i1, i3, cliff3, true);
                }
              } else {
                // clifford product
                double cliff1 = omega_table_4[i3 + nHL * (i2 + nHL * (i1 + nHL * k))].real();
                double cliff2 = omega_table_4[i2 + nHL * (i3 + nHL * (i1 + nHL * k))].real();
                double cliff3 = omega_table_4[i3 + nHL * (i1 + nHL * (i2 + nHL * k))].real();

                if (fabs(cliff1) > 1e-10) {
                  res += compute_B4(k, i1, i2, i3, cliff1, false);
                  res += compute_B4(k, i1, i3, i2, cliff2, false);
                  res += compute_B4(k, i2, i1, i3, cliff3, false);
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
  for (int i = 0; i < nHL; ++i) {
    if (i != k) {
      res += compute_B2(k, i);
    }
  }

  // all indices equal
  res += compute_B(k);

  if (herm) {
    return 2 * (res + res.t());
  } else {
    return 4 * res;
  }
}

cx_mat DiracOperator::compute_B4(const int& k,
                                 const int& i2,
                                 const int& i3,
                                 const int& i4,
                                 const double& cliff,
                                 const bool& neg) const {
  // base matrix products
  cx_mat M2M3 = mat[i2] * mat[i3];
  cx_mat M2M4 = mat[i2] * mat[i4];
  cx_mat M3M4 = mat[i3] * mat[i4];
  cx_mat M2M3M4 = M2M3 * mat[i4];

  // return value
  cx_mat res(dim, dim, fill::eye);

  if (neg) {

    // traces
    double tr234 = trace(M2M3M4).imag();
    double tr2 = trace(mat[i2]).real();
    double tr3 = trace(mat[i3]).real();
    double tr4 = trace(mat[i4]).real();

    // compute sum
    res *= -2 * eps[k] * tr234;
    res += cx_double(0., dim) * (M2M3M4 - M2M3M4.t());
    res += cx_double(0., eps[i2] * tr2) * (M3M4 - M3M4.t());
    res += cx_double(0., eps[i3] * tr3) * (M2M4 - M2M4.t());
    res += cx_double(0., eps[i4] * tr4) * (M2M3 - M2M3.t());
  } else {
    // traces
    double tr234 = trace(M2M3M4).real();
    double tr23 = trace(M2M3).real();
    double tr24 = trace(M2M4).real();
    double tr34 = trace(M3M4).real();
    double tr2 = trace(mat[i2]).real();
    double tr3 = trace(mat[i3]).real();
    double tr4 = trace(mat[i4]).real();

    // compute sum
    res *= 2 * eps[k] * tr234;
    res += dim * (M2M3M4 + M2M3M4.t());
    res += eps[i2] * tr2 * (M3M4 + M3M4.t());
    res += eps[i3] * tr3 * (M2M4 + M2M4.t());
    res += eps[i4] * tr4 * (M2M3 + M2M3.t());
    res += 2 * eps[k] * eps[i2] * tr34 * mat[i2];
    res += 2 * eps[k] * eps[i3] * tr24 * mat[i3];
    res += 2 * eps[k] * eps[i4] * tr23 * mat[i4];
  }

  return cliff * res;
}

cx_mat DiracOperator::compute_B2(const int& k, const int& i) const {
  // clifford product
  double cliff = omega_table_4[i + nHL * (k + nHL * (i + nHL * k))].real();

  // base matrix products
  cx_mat MiMk = mat[i] * mat[k];
  cx_mat MiMi = mat[i] * mat[i];
  cx_mat MiMiMk = mat[i] * MiMk;
  cx_mat MiMkMi = MiMk * mat[i];

  // traces
  double triki = trace(MiMkMi).real();
  double trik = trace(MiMk).real();
  double trii = trace(MiMi).real();
  double tri = trace(mat[i]).real();
  double trk = trace(mat[k]).real();

  // return value
  cx_mat res(dim, dim, fill::eye);

  if (cliff < 0) {
    // compute sum
    res *= eps[k] * triki;
    res += dim * (MiMiMk + MiMiMk.t() - MiMkMi);
    res += eps[i] * tri * (MiMk + MiMk.t());
    res += 2 * eps[k] * eps[i] * trik * mat[i];
    res += eps[k] * trk * MiMi;
    res += trii * mat[k];
  } else {
    // compute sum
    res *= 3 * eps[k] * triki;
    res += dim * (MiMiMk + MiMiMk.t() + MiMkMi);
    res += 3 * eps[i] * tri * (MiMk + MiMk.t());
    res += 6 * eps[k] * eps[i] * trik * mat[i];
    res += 3 * eps[k] * trk * MiMi;
    res += 3 * trii * mat[k];
  }

  return 2 * dim_omega * res;
}

cx_mat DiracOperator::compute_B(const int& k) const {
  // base matrix products
  cx_mat M2 = mat[k] * mat[k];
  cx_mat M3 = mat[k] * M2;

  // traces
  double tr3 = trace(M3).real();
  double tr2 = trace(M2).real();
  double tr1 = trace(mat[k]).real();

  cx_mat res(dim, dim, fill::eye);
  res *= eps[k] * tr3;
  res += dim * M3;
  res += 3 * tr2 * mat[k];
  res += 3 * eps[k] * tr1 * M2;

  return 2 * dim_omega * res;
}

void DiracOperator::randomise(gsl_rng* engine) {
  for (int i = 0; i < nHL; ++i) {
    // loop on indices
    for (int j = 0; j < dim; ++j) {
      double x;
      x = gsl_ran_gaussian(engine, 1.);
      mat[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < dim; ++k) {
        double a, b;
        a = gsl_ran_gaussian(engine, 1.);
        b = gsl_ran_gaussian(engine, 1.);
        mat[i](j, k) = cx_double(a, b) / sqrt(2.);
        mat[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }

  for (int i = 0; i < nHL; ++i) {
    // loop on indices
    for (int j = 0; j < dim; ++j) {
      double x;
      x = gsl_ran_gaussian(engine, 1.);
      mom[i](j, j) = cx_double(x, 0.);

      for (int k = j + 1; k < dim; ++k) {
        double a, b;
        a = gsl_ran_gaussian(engine, 1.);
        b = gsl_ran_gaussian(engine, 1.);
        mom[i](j, k) = cx_double(a, b) / sqrt(2.);
        mom[i](k, j) = cx_double(a, -b) / sqrt(2.);
      }
    }
  }
}
