//
// Created by Paul Druce on 12/11/2022.
//
#include "Action.hpp"

using namespace std;
using namespace arma;

Action::Action(double g2_, double g4_)
	: g2(g2_), g4(g4_) {}

double Action::calculate_S_from_dirac(const DiracOperator &D) const {
  cx_mat dirac = D.build_dirac();
  cx_mat dirac2 = dirac * dirac;
  double trdirac2 = trace(dirac2).real();
  double trdirac4 = trace(dirac2 * dirac2).real();
  return g2 * trdirac2 + g4 * trdirac4;
}

double Action::calculate_S(const DiracOperator &D) const {
  return g2 * dirac2(D) + g4 * dirac4(D);
}

double Action::dirac2(const DiracOperator &D) const {
  auto *mat = D.get_mats();
  auto *eps = D.get_eps();
  double res = 0.;
  for (int i = 0; i < D.nHL; ++i) {
	double tr2 = trace(mat[i] * mat[i]).real();
	double tr1 = trace(mat[i]).real();

	res += (D.dim * tr2 + eps[i] * tr1 * tr1);
  }

  return 2. * D.dim_omega * res;
}

double Action::dirac4(const DiracOperator &D) const {
  double res = 0.;

  // four distinct indices
  for (int i = 0; i < D.nHL; ++i) {
	for (int j = i + 1; j < D.nHL; ++j) {
	  for (int k = j + 1; k < D.nHL; ++k) {
		for (int l = k + 1; l < D.nHL; ++l)
		  res += 8 * (compute_A4(D, i, j, k, l) + compute_A4(D, i, j, l, k) + compute_A4(D, i, k, j, l));
	  }
	}
  }

  // two distinct pairs of equal indices
  for (int i = 0; i < D.nHL; ++i) {
	for (int j = i + 1; j < D.nHL; ++j)
	  res += 2 * compute_A2(D, i, j);
  }

  // all indices equal
  for (int i = 0; i < D.nHL; ++i)
	res += compute_A(D, i);

  return res;
}

double Action::compute_A4(const DiracOperator &D, const int &i1, const int &i2, const int &i3, const int &i4) const {
  // epsilon factor
  auto *eps = D.get_eps();
  auto *omega_table_4 = D.get_omega_table_4();
  auto *mat = D.get_mats();

  int e = eps[i1] * eps[i2] * eps[i3] * eps[i4];

  // if e=-1, then [1+*e] becomes 2i*imag
  // and the clifford part is guaranteed to
  // be pure imaginary
  if (e < 0) {
	// clifford product
	double cliff = omega_table_4[i4 + D.nHL * (i3 + D.nHL * (i2 + D.nHL * i1))].imag();

	if (fabs(cliff) > 1e-10) {
	  // base matrix products
	  cx_mat M1M2 = mat[i1] * mat[i2];
	  cx_mat M1M3 = mat[i1] * mat[i3];
	  cx_mat M1M4 = mat[i1] * mat[i4];
	  cx_mat M2M3 = mat[i2] * mat[i3];
	  cx_mat M2M4 = mat[i2] * mat[i4];
	  cx_mat M3M4 = mat[i3] * mat[i4];

	  // traces
	  double tr1234 = trace(M1M2 * M3M4).imag();
	  double tr234 = trace(M2M3 * mat[i4]).imag();
	  double tr134 = trace(M1M3 * mat[i4]).imag();
	  double tr124 = trace(M1M2 * mat[i4]).imag();
	  double tr123 = trace(M1M2 * mat[i3]).imag();
	  double tr1 = trace(mat[i1]).real();
	  double tr2 = trace(mat[i2]).real();
	  double tr3 = trace(mat[i3]).real();
	  double tr4 = trace(mat[i4]).real();

	  // compute sum
	  double res = D.dim * tr1234;
	  res += eps[i1] * tr1 * tr234;
	  res += eps[i2] * tr2 * tr134;
	  res += eps[i3] * tr3 * tr124;
	  res += eps[i4] * tr4 * tr123;

	  return -2 * cliff * res;
	  // NOTE: this minus here comes from the 'i' in cliff
	  // and the 'i' coming from 2i*imag
	} else {
	  return 0.;
	}
  } else {
	// clifford product
	double cliff = omega_table_4[i4 + D.nHL * (i3 + D.nHL * (i2 + D.nHL * i1))].real();

	if (fabs(cliff) > 1e-10) {
	  // base matrix products
	  cx_mat M1M2 = mat[i1] * mat[i2];
	  cx_mat M1M3 = mat[i1] * mat[i3];
	  cx_mat M1M4 = mat[i1] * mat[i4];
	  cx_mat M2M3 = mat[i2] * mat[i3];
	  cx_mat M2M4 = mat[i2] * mat[i4];
	  cx_mat M3M4 = mat[i3] * mat[i4];

	  // traces
	  double tr1234 = trace(M1M2 * M3M4).real();
	  double tr234 = trace(M2M3 * mat[i4]).real();
	  double tr134 = trace(M1M3 * mat[i4]).real();
	  double tr124 = trace(M1M2 * mat[i4]).real();
	  double tr123 = trace(M1M2 * mat[i3]).real();
	  double tr12 = trace(M1M2).real();
	  double tr34 = trace(M3M4).real();
	  double tr13 = trace(M1M3).real();
	  double tr24 = trace(M2M4).real();
	  double tr14 = trace(M1M4).real();
	  double tr23 = trace(M2M3).real();
	  double tr1 = trace(mat[i1]).real();
	  double tr2 = trace(mat[i2]).real();
	  double tr3 = trace(mat[i3]).real();
	  double tr4 = trace(mat[i4]).real();

	  double res = D.dim * tr1234;
	  res += eps[i1] * tr1 * tr234;
	  res += eps[i2] * tr2 * tr134;
	  res += eps[i3] * tr3 * tr124;
	  res += eps[i4] * tr4 * tr123;
	  res += eps[i1] * eps[i2] * tr12 * tr34;
	  res += eps[i1] * eps[i3] * tr13 * tr24;
	  res += eps[i1] * eps[i4] * tr14 * tr23;

	  cliff = omega_table_4[i4 + D.nHL * (i3 + D.nHL * (i2 + D.nHL * i1))].real();

	  return 2 * cliff * res;
	} else {
	  return 0.;
	}
  }
}

double Action::compute_A2(const DiracOperator &D, const int &i1, const int &i2) const {
  auto *omega_table_4 = D.get_omega_table_4();
  auto *mat = D.get_mats();
  auto *eps = D.get_eps();

  // clifford product
  double cliff = omega_table_4[i2 + D.nHL * (i1 + D.nHL * (i2 + D.nHL * i1))].real();

  // base matrix products
  cx_mat M1M1 = mat[i1] * mat[i1];
  cx_mat M2M2 = mat[i2] * mat[i2];
  cx_mat M1M2 = mat[i1] * mat[i2];

  // traces
  double tr1122 = trace(M1M1 * M2M2).real();
  double tr1212 = trace(M1M2 * M1M2).real();
  double tr122 = trace(M1M2 * mat[i2]).real();
  double tr112 = trace(M1M1 * mat[i2]).real();
  double tr11 = trace(M1M1).real();
  double tr22 = trace(M2M2).real();
  double tr12 = trace(M1M2).real();
  double tr1 = trace(mat[i1]).real();
  double tr2 = trace(mat[i2]).real();

  if (cliff < 0) {
	// compute sum
	double res = D.dim * (2 * tr1122 - tr1212);
	res += 2 * eps[i1] * tr1 * tr122;
	res += 2 * eps[i2] * tr2 * tr112;
	res += tr11 * tr22;
	res += 2 * eps[i1] * eps[i2] * tr12 * tr12;

	return 2 * D.dim_omega * res;
  } else {
	// compute sum
	double res = D.dim * (2 * tr1122 + tr1212);
	res += 6 * eps[i1] * tr1 * tr122;
	res += 6 * eps[i2] * tr2 * tr112;
	res += 3 * tr11 * tr22;
	res += 6 * eps[i1] * eps[i2] * tr12 * tr12;

	return 2 * D.dim_omega * res;
  }
}

double Action::compute_A(const DiracOperator &D, const int &i) const {
  auto *mat = D.get_mats();
  auto *eps = D.get_eps();

  // base matrix products
  cx_mat M2 = mat[i] * mat[i];
  cx_mat M3 = M2 * mat[i];

  // traces
  double tr1 = trace(mat[i]).real();
  double tr2 = trace(M2).real();
  double tr3 = trace(M3).real();
  double tr4 = trace(M3 * mat[i]).real();

  double res = D.dim * tr4;
  res += 4 * eps[i] * tr1 * tr3;
  res += 3 * tr2 * tr2;

  return 2 * D.dim_omega * res;
}

void Action::set_params(double g2_, double g4_) {
  this->g2 = g2_;
  this->g4 = g4_;
}
void Action::set_g4(double value) { this->g4 = value; }
void Action::set_g2(double value) { this->g2 = value; }



