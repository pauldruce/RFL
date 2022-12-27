//
// Created by Paul Druce on 13/11/2022.
//
#include "Metropolis.hpp"
#include <cmath>

using namespace std;
using namespace arma;

double Metropolis::delta24(const DiracOperator &D,
						   const Action &A,
						   const int &x,
						   const int &I,
						   const int &J,
						   const cx_double &z) const {
  return A.get_g2() * delta2(D, A, x, I, J, z) + delta4(D, A, x, I, J, z);
}

double Metropolis::delta2(const DiracOperator &D,
						  const Action &A,
						  const int &x,
						  const int &I,
						  const int &J,
						  const cx_double &z) const {
  auto *mat = D.get_mats();
  auto *eps = D.get_eps();

  if (I != J) {
	return 4. * D.dim_omega * D.dim * (2. * (z * mat[x](J, I)).real() + norm(z));
  } else {
	double trM = trace(mat[x]).real();
	return 8. * D.dim_omega * z.real() * (D.dim * (mat[x](I, I).real() + z.real()) + eps[x] * (trM + z.real()));
  }
}

double Metropolis::delta4(const DiracOperator &D,
						  const Action &A,
						  const int &x,
						  const int &I,
						  const int &J,
						  const cx_double &z) const {
  double res = 0.;

  auto *omega_table_4 = D.get_omega_table_4();
  auto *mat = D.get_mats();
  auto *eps = D.get_eps();

  // D^3 dD part
  for (int i3 = 0; i3 < D.nHL; ++i3) {
	for (int i2 = 0; i2 < D.nHL; ++i2) {
	  for (int i1 = 0; i1 <= i3; ++i1) {
		cx_double cliff = omega_table_4[x + D.nHL * (i3 + D.nHL * (i2 + D.nHL * i1))];

		if (fabs(cliff.real()) > 1e-10 || fabs(cliff.imag()) > 1e-10) {
		  // compute necessary matrix products
		  cx_mat M1M2 = mat[i1] * mat[i2];
		  cx_mat M2M3 = mat[i2] * mat[i3];
		  cx_mat M1M3 = mat[i1] * mat[i3];
		  cx_mat M1M2M3 = mat[i1] * M2M3;

		  // compute necessary traces
		  double trM1 = trace(mat[i1]).real();
		  double trM2 = trace(mat[i2]).real();
		  double trM3 = trace(mat[i3]).real();
		  double trM1M2 = trace(M1M2).real();
		  double trM2M3 = trace(M2M3).real();
		  double trM1M3 = trace(M1M3).real();
		  cx_double trM1M2M3 = trace(M1M2M3);

		  // off-diagonal update
		  if (I != J) {

			// compute terms
			// _______________________________________________________________________________________
			cx_double T1 = M1M2M3(J, I) * z + M1M2M3(I, J) * conj(z);
			T1 = T1 + conj(T1) * (double)(eps[i1] * eps[i2] * eps[i3] * eps[x]);
			T1 *= (double)D.dim;

			cx_double T2 = M1M2(J, I) * z + M1M2(I, J) * conj(z);
			T2 = T2 * (double)(eps[i3]) + conj(T2) * (double)(eps[i1] * eps[i2] * eps[x]);
			T2 = T2 * trM3;
			T1 += T2;

			cx_double T3 = M1M3(J, I) * z + M1M3(I, J) * conj(z);
			T3 = T3 * (double)(eps[i2]) + conj(T3) * (double)(eps[i1] * eps[i3] * eps[x]);
			T3 = T3 * trM2;
			T1 += T3;

			cx_double T4 = M2M3(J, I) * z + M2M3(I, J) * conj(z);
			T4 = T4 * (double)(eps[i1]) + conj(T4) * (double)(eps[i2] * eps[i3] * eps[x]);
			T4 = T4 * trM1;
			T1 += T4;

			double T5 = trM1M2 * (eps[i1] * eps[i2] + eps[i3] * eps[x]);
			T5 *= 2. * (mat[i3](J, I) * z).real();
			T1 += T5;

			double T6 = trM2M3 * (eps[i2] * eps[i3] + eps[i1] * eps[x]);
			T6 *= 2. * (mat[i1](J, I) * z).real();
			T1 += T6;

			double T7 = trM1M3 * (eps[i1] * eps[i3] + eps[i2] * eps[x]);
			T7 *= 2. * (mat[i2](J, I) * z).real();
			T1 += T7;
			//________________________________________________________________________________________



			// add to total
			if (i1 != i3)
			  res += 2. * (cliff * T1).real();
			else
			  res += (cliff * T1).real();
		  }


			// diagonal update
		  else {

			// compute terms
			// _______________________________________________________________________________________
			cx_double T1 = M1M2M3(I, I);
			T1 = T1 + conj(T1) * (double)(eps[i1] * eps[i2] * eps[i3] * eps[x]);
			T1 = T1 * (double)D.dim;

			cx_double T2 = M1M2(I, I);
			T2 = T2 * (double)(eps[i3]) + conj(T2) * (double)(eps[i1] * eps[i2] * eps[x]);
			T2 *= trM3;
			T1 += T2;

			cx_double T3 = M1M3(I, I);
			T3 = T3 * (double)(eps[i2]) + conj(T3) * (double)(eps[i1] * eps[i3] * eps[x]);
			T3 *= trM2;
			T1 += T3;

			cx_double T4 = M2M3(I, I);
			T4 = T4 * (double)(eps[i1]) + conj(T4) * (double)(eps[i2] * eps[i3] * eps[x]);
			T4 *= trM1;
			T1 += T4;

			double T5 = trM1M2 * (eps[i1] * eps[i2] + eps[i3] * eps[x]);
			T5 *= mat[i3](I, I).real();
			T1 += T5;

			double T6 = trM2M3 * (eps[i2] * eps[i3] + eps[i1] * eps[x]);
			T6 *= mat[i1](I, I).real();
			T1 += T6;

			double T7 = trM1M3 * (eps[i1] * eps[i3] + eps[i2] * eps[x]);
			T7 *= mat[i2](I, I).real();
			T1 += T7;

			cx_double T8 = conj(trM1M2M3) * (double)(eps[i1] * eps[i2] * eps[i3]) + trM1M2M3 * (double)(eps[x]);
			T1 += T8;
			//________________________________________________________________________________________



			// add to total
			if (i1 != i3)
			  res += (cliff * T1).real() * 4. * z.real();
			else
			  res += (cliff * T1).real() * 2. * z.real();
		  }
		}
	  }
	}
  }

  res *= 4.;



  // D^2 dD^2 and D dD D dD term
  double temp = 0;
  for (int i = 0; i < D.nHL; ++i) {
	double cliff = omega_table_4[x + D.nHL * (i + D.nHL * (x + D.nHL * i))].real();

	// compute necessary matrix products
	cx_mat M1M1 = mat[i] * mat[i];

	// compute necessary traces
	double trM1 = trace(mat[i]).real();
	double trM1M1 = trace(M1M1).real();

	// off-diagonal update
	if (I != J) {
	  // compute terms D^2 dD^2
	  // _______________________________________________________________________________________
	  double T11 = 2 * D.dim * (M1M1(I, I).real() + M1M1(J, J).real());
	  double T21 = 4 * eps[i] * trM1 * (mat[i](I, I).real() + mat[i](J, J).real());
	  double T31 = (z * mat[i](J, I)).real();
	  T31 *= T31 * 16 * eps[i] * eps[x];
	  //________________________________________________________________________________________

	  // compute terms D dD D dD
	  // _______________________________________________________________________________________

	  double T12 = (mat[i](J, I) * mat[i](J, I) * z * z).real();
	  T12 += mat[i](I, I).real() * mat[i](J, J).real() * norm(z);
	  T12 *= 4 * D.dim;

	  double T22 = 4 * eps[i] * trM1 * (mat[i](I, I).real() + mat[i](J, J).real());
	  double T32 = (mat[i](J, I) * z).real();
	  T32 *= T32 * 16 * eps[i] * eps[x];
	  //________________________________________________________________________________________



	  // add to total
	  temp += 2. * D.dim_omega * (norm(z) * (T11 + T21 + 4. * trM1M1) + T31);
	  temp += cliff * (T12 + norm(z) * (T22 + 4. * trM1M1) + T32);

	}

	  // diagonal update
	else {
	  // compute terms D^2 dD^2
	  // _______________________________________________________________________________________
	  double T11 = 2. * D.dim * M1M1(I, I).real();
	  double T21 = 4. * eps[x] * M1M1(I, I).real();
	  double T31 = 4. * eps[i] * trM1 * mat[i](I, I).real();
	  double T41 = mat[i](I, I).real();
	  T41 *= T41 * 4. * eps[i] * eps[x];
	  //________________________________________________________________________________________

	  // compute terms D dD D dD
	  // _______________________________________________________________________________________
	  double T12 = mat[i](I, I).real();
	  T12 *= T12 * 2. * D.dim;
	  double T22 = 4. * eps[x] * M1M1(I, I).real();
	  double T32 = 4. * eps[i] * trM1 * mat[i](I, I).real();
	  double T42 = mat[i](I, I).real();
	  T42 *= T42 * 4. * eps[i] * eps[x];
	  //________________________________________________________________________________________

	  // add to total
	  temp += 8. * z.real() * z.real() * D.dim_omega * (T11 + T21 + T31 + T41 + 2. * trM1M1);
	  temp += 4. * z.real() * z.real() * cliff * (T12 + T22 + T32 + T42 + 2. * trM1M1);
	}
  }

  res += 2. * temp;



  // D dD^3 term

  // off-diagonal update
  if (I != J) {
	temp = 4. * D.dim_omega * (D.dim + 6) * norm(z) * (mat[x](J, I) * z).real();
	res += 4. * temp;
  }

	// diagonal update
  else {
	double trMx = trace(mat[x]).real();
	double rez = 2. * z.real();
	temp = 2. * rez * rez * rez * D.dim_omega * (mat[x](I, I).real() * (D.dim + 3. * eps[x] + 3.) + eps[x] * trMx);
	res += 4. * temp;
  }


  // dD^4 term

  // off-diagonal update
  if (I != J) {
	temp = D.dim_omega * 4. * norm(z) * norm(z) * (D.dim + 6.);
	res += temp;
  }

	// diagonal update
  else {
	double rez = z.real();
	temp = D.dim_omega * 32. * (D.dim + 3. + 4 * eps[x]) * rez * rez * rez * rez;
	res += temp;
  }

  return res;
}

void Metropolis::MMC_duav(const DiracOperator &D,
						  const Action &A,
						  double &scale,
						  const int &iter,
						  gsl_rng *engine,
						  const double &target) const {
  // initial (_i) and final (_f) action2 and action4
  auto *s_i = new double[2];
  auto *s_f = new double[2];

  // calculate length of a sweep in terms of dofs
  int Nsw = D.nHL * D.dim * D.dim - D.nL;

  // dual averaging variables
  const double shr = 0.05;
  const double kappa = 0.75;
  const int i0 = 10;
  double Stat = 0;
  double mu = log(10 * scale);
  double log_scale_avg = log(scale);

  // iter sweeps of metropolis
  for (int i = 0; i < iter; ++i) {
	for (int j = 0; j < Nsw; ++j) {
	  // set action to previous final value,
	  // unless it's the first iteration
	  if (j) {
		s_i[0] = s_f[0];
		s_i[1] = s_f[1];
	  } else {
		s_i[0] = A.dirac2(D);
		s_i[1] = A.dirac4(D);
	  }

	  Stat += target - MMC_duav_core(D, A, scale, engine, s_i, s_f);

	  // perform dual averaging
	  double log_scale = mu - Stat * sqrt(i + 1) / (shr * (i + 1 + i0));
	  scale = exp(log_scale);
	  double eta = pow(i + 1, -kappa);
	  log_scale_avg = eta * log_scale + (1 - eta) * log_scale_avg;
	}
  }

  // set scale on its final dual averaged value
  scale = exp(log_scale_avg);

  delete[] s_i;
  delete[] s_f;
}

double Metropolis::MMC(const DiracOperator &D,
					   const Action &A,
					   const double &scale,
					   const int &iter,
					   gsl_rng *engine) const {
  // initial (_i) and final (_f) action2 and action4
  auto *s_i = new double[2];
  auto *s_f = new double[2];

  // calculate length of a sweep in terms of dofs
  int Nsw = D.nHL * D.dim * D.dim - D.nL;

  // return statistic
  double Stat = 0;

  // iter sweeps of metropolis
  for (int i = 0; i < iter; ++i) {
	for (int j = 0; j < Nsw; ++j) {
	  // set action to previous final value,
	  // unless it's the first iteration
	  if (j) {
		s_i[0] = s_f[0];
		s_i[1] = s_f[1];
	  } else {
		s_i[0] = A.dirac2(D);
		s_i[1] = A.dirac4(D);
	  }

	  Stat += MMC_core(D, A, scale, engine, s_i, s_f);
	}
  }

  delete[] s_i;
  delete[] s_f;

  return (Stat / (iter * Nsw));
}

double Metropolis::MMC_duav_core(const DiracOperator &D,
								 const Action &A,
								 const double &scale,
								 gsl_rng *engine,
								 double *s_i,
								 double *s_f) const {
  // acceptance probability
  double e;

  // metropolis
  int x = D.nHL * gsl_rng_uniform(engine);
  int I = D.dim * gsl_rng_uniform(engine);
  int J = D.dim * gsl_rng_uniform(engine);
  double re = 0;
  double im = 0;
  cx_double z;
  if (I != J) {
	re = scale * (-1. + 2. * gsl_rng_uniform(engine));
	im = scale * (-1. + 2. * gsl_rng_uniform(engine));
	z = cx_double(re, im);
  } else {
	re = scale * (-1. + 2. * gsl_rng_uniform(engine));
	z = cx_double(re, 0);
  }

  double dS2 = delta2(D, A, x, I, J, z);
  double dS4 = delta4(D, A, x, I, J, z);
  double dS = A.get_g2() * dS2 + dS4;

  auto *mat = D.get_mats();
  // metropolis test
  if (dS < 0) {
	// update matrix element
	if (I != J) {
	  mat[x](I, J) += z;
	  mat[x](J, I) += conj(z);
	} else
	  mat[x](I, I) += 2. * z;

	// update action
	s_f[0] = s_i[0] + dS2;
	s_f[1] = s_i[1] + dS4;

	// move accepted
	e = 1;
  } else {
	e = exp(-dS);
	double p = gsl_rng_uniform(engine);

	if (e > p) {
	  // update matrix element
	  if (I != J) {
		mat[x](I, J) += z;
		mat[x](J, I) += conj(z);
	  } else
		mat[x](I, I) += 2. * z;

	  // update action
	  s_f[0] = s_i[0] + dS2;
	  s_f[1] = s_i[1] + dS4;
	} else {
	  s_f[0] = s_i[0];
	  s_f[1] = s_i[1];
	}
  }

  return e;
}

double Metropolis::MMC_core(const DiracOperator &D,
							const Action &A,
							const double &scale,
							gsl_rng *engine,
							double *s_i,
							double *s_f) const {
  // acceptance probability
  double ret = 0;

  // metropolis
  int x = D.nHL * gsl_rng_uniform(engine);
  int I = D.dim * gsl_rng_uniform(engine);
  int J = D.dim * gsl_rng_uniform(engine);
  double re = 0;
  double im = 0;
  cx_double z;
  if (I != J) {
	re = scale * (-1. + 2. * gsl_rng_uniform(engine));
	im = scale * (-1. + 2. * gsl_rng_uniform(engine));
	z = cx_double(re, im);
  } else {
	re = scale * (-1. + 2. * gsl_rng_uniform(engine));
	z = cx_double(re, 0);
  }

  double dS2 = delta2(D, A, x, I, J, z);
  double dS4 = delta4(D, A, x, I, J, z);
  double dS = A.get_g2() * dS2 + dS4;

  auto *mat = D.get_mats();
  // metropolis test
  if (dS < 0) {
	// update matrix element
	if (I != J) {
	  mat[x](I, J) += z;
	  mat[x](J, I) += conj(z);
	} else
	  mat[x](I, I) += 2. * z;

	// update action
	s_f[0] = s_i[0] + dS2;
	s_f[1] = s_i[1] + dS4;

	// move accepted
	ret = 1;
  } else {
	double e = exp(-dS);
	double p = gsl_rng_uniform(engine);

	if (e > p) {
	  // update matrix element
	  if (I != J) {
		mat[x](I, J) += z;
		mat[x](J, I) += conj(z);
	  } else
		mat[x](I, I) += 2. * z;

	  // update action
	  s_f[0] = s_i[0] + dS2;
	  s_f[1] = s_i[1] + dS4;

	  // move accepted
	  ret = 1;
	} else {
	  s_f[0] = s_i[0];
	  s_f[1] = s_i[1];
	}
  }

  return ret;
}

