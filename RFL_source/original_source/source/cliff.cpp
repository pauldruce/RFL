#include "Cliff.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace arma;

// Constructors

Cliff::Cliff(int mode) {
  //(1,0)
  if (mode == 3) {
    p = 1;
    q = 0;
    dim_gamma = 1;

    gamma.emplace_back(1, 1, fill::eye);
    chiral = cx_mat(1, 1, fill::eye);
  }
  //(0,1)
  else if (mode == 4) {
    p = 0;
    q = 1;
    dim_gamma = 1;

    cx_double z(0., -1.);
    cx_mat tmp = cx_mat(1, 1);
    tmp(0, 0) = z;
    gamma.push_back(tmp);
    chiral = cx_mat(1, 1, fill::eye);
  }
  //(2,0)
  else if (mode == 0) {
    p = 2;
    q = 0;
    dim_gamma = 2;

    cx_double z(0., 1.);

    cx_mat tmp1 = cx_mat(2, 2, fill::zeros);
    tmp1(0, 0) = 1.;
    tmp1(1, 1) = -1.;
    gamma.push_back(tmp1);

    cx_mat tmp2 = cx_mat(2, 2, fill::zeros);
    tmp2(0, 1) = 1.;
    tmp2(1, 0) = 1.;
    gamma.push_back(tmp2);

    chiral = cx_mat(2, 2, fill::zeros);
    chiral(0, 1) = z;
    chiral(1, 0) = -z;
  }
  //(1,1)
  else if (mode == 2) {
    p = 1;
    q = 1;
    dim_gamma = 2;

    cx_mat tmp1 = cx_mat(2, 2, fill::zeros);
    tmp1(0, 0) = 1.;
    tmp1(1, 1) = -1.;
    gamma.push_back(tmp1);

    cx_mat tmp2 = cx_mat(2, 2, fill::zeros);
    tmp2(0, 1) = 1.;
    tmp2(1, 0) = -1.;
    gamma.push_back(tmp2);

    chiral = cx_mat(2, 2, fill::zeros);
    chiral(0, 1) = 1;
    chiral(1, 0) = 1;
  }
  //(0,2)
  else if (mode == 1) {
    p = 0;
    q = 2;
    dim_gamma = 2;

    cx_double z(0., 1.);

    cx_mat tmp1 = cx_mat(2, 2, fill::zeros);
    tmp1(0, 0) = z;
    tmp1(1, 1) = -z;
    gamma.push_back(tmp1);

    cx_mat tmp2 = cx_mat(2, 2, fill::zeros);
    tmp2(0, 1) = 1.;
    tmp2(1, 0) = -1.;
    gamma.push_back(tmp2);

    chiral = cx_mat(2, 2, fill::zeros);
    chiral(0, 1) = 1.;
    chiral(1, 0) = 1.;
  }
}

Cliff::Cliff(int p_, int q_)
    : p(p_), q(q_), dim_gamma() {
  //(1,0)
  if (p == 1 && q == 0) {
    (*this) = Cliff(3);
  }
  //(0,1)
  if (p == 0 && q == 1) {
    (*this) = Cliff(4);
  }
  //(2,0)
  if (p == 2 && q == 0) {
    (*this) = Cliff(0);
  }
  //(1,1)
  if (p == 1 && q == 1) {
    (*this) = Cliff(2);
  }
  //(0,2)
  if (p == 0 && q == 2) {
    (*this) = Cliff(1);
    //any other case
  } else {
    init_gamma();
    sort_gamma();
  }
}

// Copy constructor
Cliff::Cliff(const Cliff& C) {
  // copy parameters
  p = C.get_p();
  q = C.get_q();
  dim_gamma = C.get_dim_gamma();

  // copy matrices
  for (int i = 0; i < p + q; i++)
    gamma.push_back(C.get_gamma(i));

  chiral = C.get_chiral();
}

// Operator =
Cliff& Cliff::operator=(const Cliff& C) {
  p = C.get_p();
  q = C.get_q();
  dim_gamma = C.get_dim_gamma();

  // delete, reallocate and copy matrices
  gamma.clear();
  for (int i = 0; i < p + q; i++)
    gamma.push_back(C.get_gamma(i));

  chiral = C.get_chiral();

  return *this;
}

static void decomp(const int p, const int q, int* dec) {
  if (p) {
    if (!(p % 2)) {
      dec[0] = p / 2;
      dec[3] = 0;
    } else {
      dec[0] = (p - 1) / 2;
      dec[3] = 1;
    }
  } else {
    dec[0] = 0;
    dec[3] = 0;
  }

  if (q) {
    if (!(q % 2)) {
      dec[1] = q / 2;
      dec[4] = 0;
    } else {
      dec[1] = (q - 1) / 2;
      dec[4] = 1;
    }
  } else {
    dec[1] = 0;
    dec[4] = 0;
  }

  if (dec[3] && dec[4]) {
    dec[3] = 0;
    dec[4] = 0;
    dec[2] = 1;
  } else {
    dec[2] = 0;
  }
}

// init_gamma gets called only if p+q > 2
void Cliff::init_gamma() {
  int dec[5];
  decomp(p, q, dec);

  vector<Cliff> vec;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < dec[i]; ++j)
      vec.emplace_back(i);
  }

  const auto begin = vec.begin();
  const auto end = vec.end();

  Cliff C1 = *begin;

  for (auto iter = begin + 1; iter != end; ++iter)
    C1 *= *iter;

  *this = C1;
}

Cliff& Cliff::operator*=(const Cliff& C2) {
  // store C2 frequently used variables
  const int p2 = C2.get_p();
  const int q2 = C2.get_q();
  const int dim2 = C2.get_dim_gamma();

  // temporary variables to avoid overwriting on (*this)
  vector<cx_mat> gamma_;
  const int p_ = p + p2;
  const int q_ = q + q2;
  const int dim_gamma_ = dim_gamma * dim2;

  // start computing product
  const cx_mat id2(dim2, dim2, fill::eye);

  gamma_.reserve(p + q);
  for (int i = 0; i < p + q; ++i)
    gamma_.emplace_back(kron(gamma[i], id2));
  for (int i = 0; i < p2 + q2; ++i)
    gamma_.emplace_back(kron(chiral, C2.get_gamma(i)));

  // compute chirality
  const int s2 = (q2 - p2 + 8 * p2) % 8;// +8*p2 is necessary becase % does not mean modulo for negative numbers
  if ((s2 % 8) % 2) {
    const int s = (q - p + 8 * p) % 8;
    const cx_mat id1(dim_gamma, dim_gamma, fill::eye);
    if (s == 2 || s == 6) {
      chiral = kron(-1 * id1, C2.get_chiral());
    } else {
      chiral = kron(id1, C2.get_chiral());
    }
  } else {
    chiral = kron(chiral, C2.get_chiral());
  }

  // overwrite on (*this)
  p = p_;
  q = q_;
  dim_gamma = dim_gamma_;
  gamma.clear();
  //   for (auto iter = gamma_.begin(); iter != gamma_.end(); ++iter)
  //      gamma.push_back((*iter));
  for (const auto& v : gamma_)
    gamma.push_back(v);

  return *this;
}

ostream& operator<<(ostream& out, const Cliff& C) {
  out << "Clifford (p, q) = (" << C.get_p() << ", " << C.get_q() << ") ";

  return out;
}

static bool areHermitian(const cx_mat& M1, const cx_mat& M2) {
  return (!(M2.is_hermitian()) && M1.is_hermitian());
}

void Cliff::sort_gamma() {
  sort(gamma.begin(), gamma.end(), areHermitian);
}
