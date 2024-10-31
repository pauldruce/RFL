#include "Clifford.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace arma;

// Constructors

Clifford::Clifford(int mode) {
  //(1,0)
  if (mode == 3) {
    m_p = 1;
    m_q = 0;
    m_dim_gamma = 1;

    m_gammas.emplace_back(1, 1, fill::eye);
    m_chiral = cx_mat(1, 1, fill::eye);
  }
  //(0,1)
  else if (mode == 4) {
    m_p = 0;
    m_q = 1;
    m_dim_gamma = 1;

    cx_double z(0., -1.);
    cx_mat tmp = cx_mat(1, 1);
    tmp(0, 0) = z;
    m_gammas.push_back(tmp);
    m_chiral = cx_mat(1, 1, fill::eye);
  }
  //(2,0)
  else if (mode == 0) {
    m_p = 2;
    m_q = 0;
    m_dim_gamma = 2;

    cx_double z(0., 1.);

    cx_mat tmp_1 = cx_mat(2, 2, fill::zeros);
    tmp_1(0, 0) = 1.;
    tmp_1(1, 1) = -1.;
    m_gammas.push_back(tmp_1);

    cx_mat tmp_2 = cx_mat(2, 2, fill::zeros);
    tmp_2(0, 1) = 1.;
    tmp_2(1, 0) = 1.;
    m_gammas.push_back(tmp_2);

    m_chiral = cx_mat(2, 2, fill::zeros);
    m_chiral(0, 1) = z;
    m_chiral(1, 0) = -z;
  }
  //(1,1)
  else if (mode == 2) {
    m_p = 1;
    m_q = 1;
    m_dim_gamma = 2;

    cx_mat tmp_1 = cx_mat(2, 2, fill::zeros);
    tmp_1(0, 0) = 1.;
    tmp_1(1, 1) = -1.;
    m_gammas.push_back(tmp_1);

    cx_mat tmp_2 = cx_mat(2, 2, fill::zeros);
    tmp_2(0, 1) = 1.;
    tmp_2(1, 0) = -1.;
    m_gammas.push_back(tmp_2);

    m_chiral = cx_mat(2, 2, fill::zeros);
    m_chiral(0, 1) = 1;
    m_chiral(1, 0) = 1;
  }
  //(0,2)
  else if (mode == 1) {
    m_p = 0;
    m_q = 2;
    m_dim_gamma = 2;

    cx_double z(0., 1.);

    cx_mat tmp_1 = cx_mat(2, 2, fill::zeros);
    tmp_1(0, 0) = z;
    tmp_1(1, 1) = -z;
    m_gammas.push_back(tmp_1);

    cx_mat tmp_2 = cx_mat(2, 2, fill::zeros);
    tmp_2(0, 1) = 1.;
    tmp_2(1, 0) = -1.;
    m_gammas.push_back(tmp_2);

    m_chiral = cx_mat(2, 2, fill::zeros);
    m_chiral(0, 1) = 1.;
    m_chiral(1, 0) = 1.;
  } else {
    std::string error_message = "Invalid Clifford mode entered: mode = " + std::to_string(mode);
    throw std::runtime_error(error_message);
  }
}

Clifford::Clifford(int p, int q)
    : m_p(p), m_q(q), m_dim_gamma(0) {
  //(1,0)
  if (m_p == 1 && m_q == 0) {
    *this = Clifford(3);
  }
  //(0,1)
  if (m_p == 0 && m_q == 1) {
    *this = Clifford(4);
  }
  //(2,0)
  if (m_p == 2 && m_q == 0) {
    *this = Clifford(0);
  }
  //(1,1)
  if (m_p == 1 && m_q == 1) {
    *this = Clifford(2);
  }
  //(0,2)
  if (m_p == 0 && m_q == 2) {
    *this = Clifford(1);
    //any other case
  } else {
    initGammas();
    sortGammas();
  }
}

// Copy constructor
Clifford::Clifford(const Clifford& clifford_to_copy) {
  // copy parameters
  m_p = clifford_to_copy.getP();
  m_q = clifford_to_copy.getQ();
  m_dim_gamma = clifford_to_copy.getGammaDimension();

  // copy matrices
  for (int i = 0; i < m_p + m_q; i++)
    m_gammas.push_back(clifford_to_copy.getGammaAtIndex(i));

  m_chiral = clifford_to_copy.getChiral();
}

// Operator =
Clifford& Clifford::operator=(const Clifford& clifford_to_copy) {
  m_p = clifford_to_copy.getP();
  m_q = clifford_to_copy.getQ();
  m_dim_gamma = clifford_to_copy.getGammaDimension();

  // delete, reallocate and copy matrices
  m_gammas.clear();
  for (int i = 0; i < m_p + m_q; i++)
    m_gammas.push_back(clifford_to_copy.getGammaAtIndex(i));

  m_chiral = clifford_to_copy.getChiral();

  return *this;
}

/**
 * A method to split (p,q) into a combination of the base modes.
 *
 * e.g. (5,3) = (2,0) + (2,0) + (0,2) + (1,1)
 * or   (7,2) = (2,0) + (2,0) + (2,0) + (1,0)
 * @param p
 * @param q
 * @param dec
 */
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
void Clifford::initGammas() {

  // Decompose the (p,q) into products of the base 5 types (1,0), (0,1),
  // (2,0), (1,1) and (0,2)
  int dec[5];
  decomp(m_p, m_q, dec);

  vector<Clifford> vec;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < dec[i]; ++j) {
      vec.emplace_back(i);
    }
  }

  const auto begin = vec.begin();
  const auto end = vec.end();
  Clifford c_1 = *begin;
  for (auto iter = begin + 1; iter != end; ++iter) {
    c_1 *= *iter;
  }
  *this = c_1;
}

// TODO: Go through this and make sure it agrees with Lawson+Michelson
Clifford& Clifford::operator*=(const Clifford& clifford_2) {
  // store C2 frequently used variables
  const int p_2 = clifford_2.getP();
  const int q_2 = clifford_2.getQ();
  const int dim_2 = clifford_2.getGammaDimension();

  // temporary variables to avoid overwriting on (*this)
  vector<cx_mat> gamma;

  const int p = m_p + p_2;
  const int q = m_q + q_2;
  const int dim_gamma = m_dim_gamma * dim_2;

  // start computing product
  const cx_mat id_2(dim_2, dim_2, fill::eye);

  gamma.reserve(m_p + m_q);
  for (int i = 0; i < m_p + m_q; ++i)
    gamma.emplace_back(kron(m_gammas[i], id_2));
  for (int i = 0; i < p_2 + q_2; ++i)
    gamma.emplace_back(kron(m_chiral, clifford_2.getGammaAtIndex(i)));

  // compute chirality
  const int s_2 = (q_2 - p_2 + 8 * p_2) % 8;// +8*p2 is necessary becase % does not mean modulo for negative numbers
  if (const bool s_2_even = s_2 % 8 % 2 == 0; !s_2_even) {
    const int s = (m_q - m_p + 8 * m_p) % 8;
    const cx_mat id_1(m_dim_gamma, m_dim_gamma, fill::eye);
    if (s == 2 || s == 6) {
      m_chiral = kron(-1 * id_1, clifford_2.getChiral());
    } else {
      m_chiral = kron(id_1, clifford_2.getChiral());
    }
  } else {
    m_chiral = kron(m_chiral, clifford_2.getChiral());
  }

  // overwrite on (*this)
  m_p = p;
  m_q = q;
  m_dim_gamma = dim_gamma;
  m_gammas.clear();

  for (const auto& v : gamma)
    m_gammas.push_back(v);

  return *this;
}

ostream& operator<<(ostream& out, const Clifford& clifford) {
  out << "Clifford (p, q) = (" << clifford.getP() << ", " << clifford.getQ() << ") ";

  return out;
}

bool areHermitian(const cx_mat& m_1, const cx_mat& m_2) {
  return !m_2.is_hermitian() && m_1.is_hermitian();
}

void Clifford::sortGammas() {
  sort(m_gammas.begin(), m_gammas.end(), areHermitian);
}
