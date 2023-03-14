#include "Cliff.hpp"
#include "Geom24.hpp"
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <vector>

// PD: Need to find homes for
// 		- shuffle
//		- copy/equals/ constructor,
// 		- &operator<<

using namespace std;
using namespace arma;

// PD : Moved into DiracOperator constructor
Geom24::Geom24(int p_, int q_, int dim_, double g2_)
    : p(p_), q(q_), dim(dim_), g2(g2_) {
  // initialize derived parameters
  derived_parameters();

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

// PD: Not moved into new lib
Geom24::Geom24(istream& in) {
  read_parameters(in);

  // initialize derived parameters
  derived_parameters();

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

// Copy constructor
Geom24::Geom24(const Geom24& G) {
  // copy parameters
  dim = G.get_dim();
  p = G.get_p();
  q = G.get_q();
  g2 = G.get_g2();
  nH = G.get_nH();
  nL = G.get_nL();
  nHL = G.get_nHL();
  dim_omega = G.get_dim_omega();

  // allocate and copy matrices
  mat = new arma::cx_mat[nHL];
  mom = new arma::cx_mat[nHL];
  eps = new int[nHL];
  omega = new arma::cx_mat[nHL];
  for (int i = 0; i < nHL; i++) {
    mat[i] = G.get_mat(i);
    mom[i] = G.get_mom(i);
    eps[i] = G.get_eps(i);
    omega[i] = G.get_omega(i);
  }

  int nHL4 = (int)pow(nHL, 4);
  omega_table_4 = new cx_double[nHL4];
  for (int i = 0; i < nHL4; i++)
    omega_table_4[i] = G.get_omega_table_4(i);
}

// Operator =
Geom24& Geom24::operator=(const Geom24& G) {
  dim = G.get_dim();
  p = G.get_p();
  q = G.get_q();
  g2 = G.get_g2();
  nH = G.get_nH();
  nL = G.get_nL();
  nHL = G.get_nHL();
  dim_omega = G.get_dim_omega();

  // delete, reallocate and copy matrices
  delete[] mat;
  delete[] mom;
  delete[] omega;
  delete[] eps;
  mat = new arma::cx_mat[nHL];
  mom = new arma::cx_mat[nHL];
  eps = new int[nHL];
  omega = new arma::cx_mat[nHL];
  for (int i = 0; i < nHL; i++) {
    mat[i] = G.get_mat(i);
    mom[i] = G.get_mat(i);
    eps[i] = G.get_eps(i);
    omega[i] = G.get_omega(i);
  }

  delete[] omega_table_4;
  int nHL4 = (int)pow(nHL, 4);
  omega_table_4 = new cx_double[nHL4];
  for (int i = 0; i < nHL4; i++)
    omega_table_4[i] = G.get_omega_table_4(i);

  return *this;
}

// Destructor
// PD: Moved into DiracOperator
Geom24::~Geom24() {
  delete[] mat;
  delete[] mom;
  delete[] eps;
  delete[] omega;
  delete[] omega_table_4;
}

// Read parameters from istream
istream& Geom24::read_parameters(istream& in) {
  if (in) {
    // read basic parameters
    in >> p >> q >> dim >> g2;

    // clear input stream state
    in.clear();
  }

  return in;
}

// PD: Moved into constructor of DiracOperator class
void Geom24::derived_parameters() {
  int n = p + q;

  // create a type (p, q) clifford module
  Cliff C(p, q);
  vector<cx_mat> gamma = C.get_gamma();
  dim_omega = C.get_dim_gamma();

  vector<cx_mat> herm;
  vector<cx_mat> anti;

  for (int i = 0; i < p; ++i)
    herm.push_back(gamma[i]);

  for (int i = 0; i < q; ++i)
    anti.emplace_back(cx_double(0, 1) * gamma[p + i]);

  int count = pow(2, n);
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
    int k = vec.size();
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

  nH = herm.size();
  nL = anti.size();
  nHL = nH + nL;

  omega = new cx_mat[nHL];
  for (int i = 0; i < nH; ++i)
    omega[i] = herm[i];
  for (int i = 0; i < nL; ++i)
    omega[nH + i] = anti[i];

  init_omega_table_4();
}

ostream& operator<<(ostream& out, const Geom24& G) {
  out << "Geometry (p, q, dim, g2) = (" << G.get_p() << ", " << G.get_q() << ", " << G.get_dim() << ", " << G.get_g2()
      << ") ";

  return out;
}

void Geom24::print_S(ostream& out) const {
  out << dirac2() << " " << dirac4() << endl;
}

void Geom24::print_HL(ostream& out) const {
  for (int i = 0; i < nHL; ++i) {
    for (int j = 0; j < dim; ++j) {
      for (int k = j; k < dim; ++k)
        out << mat[i](j, k).real() << " " << mat[i](j, k).imag() << " ";
    }
    out << endl;
  }
  out << endl;
}

void Geom24::reverse_mom() {
  for (int i = 0; i < nHL; ++i)
    mom[i] *= -1;
}

// PD - Moved into DiracOperator as method Randomise
void Geom24::shuffle(gsl_rng* engine) {
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

// PD: Not moved into new lib
istream& Geom24::read_mat(istream& in) {
  if (in) {
    // loop on matrices
    for (int i = 0; i < nHL; ++i) {
      // loop on indices
      for (int j = 0; j < dim; ++j) {
        for (int k = j; k < dim; ++k) {
          double x, y;
          in >> x >> y;

          if (j != k) {
            mat[i](j, k) = cx_double(x, y);
            mat[i](k, j) = cx_double(x, -y);
          } else {
            mat[i](j, j) = cx_double(x, 0.);
          }
        }
      }
    }

    // clear input stream state
    in.clear();
  }

  return in;
}

// PD: In DiracOperator
void Geom24::init_omega_table_4() {
  omega_table_4 = new cx_double[nHL * nHL * nHL * nHL];

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

// PD: In DiracOperator as static method
static vector<int> baseConversion(int dec, const int& base, const int& max) {
  vector<int> rem;

  while (dec) {
    rem.push_back(dec % base);
    dec /= base;
  }

  for (int i = rem.size(); i < max; i++)
    rem.push_back(0);

  reverse(rem.begin(), rem.end());

  return rem;
}

// PD: In DiracOperator
void Geom24::print_omega_table_4() const {
  const int n = pow(nHL, 4);

  for (int i = 0; i < n; ++i) {
    cx_double z = omega_table_4[i];
    if (z != cx_double(0., 0.)) {
      int e = 1;
      vector<int> prod = baseConversion(i, nHL, 4);
      //            vector<int>::const_iterator end(prod.end());
      //            for (vector<int>::const_iterator iter = prod.begin(); iter != end; ++iter) {
      //                cout << (*iter) << " ";
      //                e *= eps[(*iter)];
      //            }
      for (const auto& elem : prod) {
        cout << elem << " ";
        e *= eps[elem];
      }
      cout << " " << omega_table_4[i] << e << endl;
    }
  }
}

//void Geom24::adjust() {
//    // hermitianize
//    for (int i = 0; i < nHL; ++i)
//        mat[i] = 0.5 * (mat[i] + mat[i].t());
//
//    // tracelessitize
//    for (int i = nH; i < nHL; ++i) {
//        double tr = trace(mat[i]).real() / dim;
//        mat[i] = mat[i] - tr * cx_mat(dim, dim, fill::eye);
//    }
//}

// PD: In DiracOperator
cx_mat Geom24::build_dirac() const {
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
