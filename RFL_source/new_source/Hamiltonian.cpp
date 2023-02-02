//
// Created by Paul Druce on 12/11/2022.
//

#include "Hamiltonian.hpp"
#include <armadillo>
#include "Action.hpp"

using namespace std;
using namespace arma;

void Hamiltonian::sample_mom(const DiracOperator &D,
							 gsl_rng *engine) const {
  auto *mom = D.get_moms();
  for (int i = 0; i < D.nHL; ++i) {
	// loop on indices
	for (int j = 0; j < D.dim; ++j) {
	  double x;
	  x = gsl_ran_gaussian(engine, 1.);
	  mom[i](j, j) = cx_double(x, 0.);

	  for (int k = j + 1; k < D.dim; ++k) {
		double a, b;
		a = gsl_ran_gaussian(engine, 1.);
		b = gsl_ran_gaussian(engine, 1.);
		mom[i](j, k) = cx_double(a, b) / sqrt(2.);
		mom[i](k, j) = cx_double(a, -b) / sqrt(2.);
	  }
	}
  }
}

double Hamiltonian::calculate_K(const DiracOperator &D) const {
  double res = 0;
  auto *mom = D.get_moms();

  for (int i = 0; i < D.nHL; ++i)
	res += trace(mom[i] * mom[i]).real();

  return res / 2;
}

double Hamiltonian::calculate_H(const DiracOperator &D, const Action &A) const {
  return A.calculate_S(D) + calculate_K(D);
}

void Hamiltonian::leapfrog(const DiracOperator &D, const int &Nt, const double &dt, const double g2) const {
  auto *mat = D.get_mats();
  auto *mom = D.get_moms();

  for (int i = 0; i < D.nHL; ++i) {
	mat[i] += (dt / 2.) * mom[i];

	for (int j = 0; j < Nt - 1; ++j) {
	  mom[i] += -dt * D.der_dirac24(i, true, g2);
	  mat[i] += dt * mom[i];
	}

	mom[i] += -dt * D.der_dirac24(i, true, g2);
	mat[i] += (dt / 2.) * mom[i];
  }
}

void Hamiltonian::omelyan(const DiracOperator &D, const int &Nt, const double &dt, const double g2) const {
  double xi = 0.1931833;

  auto *mat = D.get_mats();
  auto *mom = D.get_moms();

  for (int i = 0; i < D.nHL; ++i) {
	mat[i] += xi * dt * mom[i];

	for (int j = 0; j < Nt - 1; ++j) {
	  mom[i] += -(dt / 2.) * D.der_dirac24(i, true, g2);
	  mat[i] += (1 - 2 * xi) * dt * mom[i];
	  mom[i] += -(dt / 2.) * D.der_dirac24(i, true, g2);
	  mat[i] += 2 * xi * dt * mom[i];
	}

	mom[i] += -(dt / 2.) * D.der_dirac24(i, true, g2);
	mat[i] += (1 - 2 * xi) * dt * mom[i];
	mom[i] += -(dt / 2.) * D.der_dirac24(i, true, g2);
	mat[i] += xi * dt * mom[i];
  }
}

void Hamiltonian::run_HMC_duav(const DiracOperator &D,
							   const Action &A,
							   const int &Nt,
							   double &dt,
							   const int &iter,
							   gsl_rng *engine,
							   const double &target,
							   const string &integrator) const{
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto *en_i = new double[4];
  auto *en_f = new double[4];

  // dual averaging variables for dt
  const double shr = 0.05;
  const double kappa = 0.75;
  const int i0 = 10;
  double Stat = 0;
  double mu = log(10 * dt);
  double log_dt_avg = log(dt);

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
	// if it's not the first interation set potential to
	// previous final value, otherwise compute it
	if (i) {
	  en_i[0] = en_f[0];
	  en_i[1] = en_f[1];
	} else {
	  en_i[0] = A.dirac2(D);
	  en_i[1] = A.dirac4(D);
	}


	// core part of HMC
	Stat += target - run_HMC_duav_core(D, A, Nt, dt, engine, en_i, en_f, integrator);

	// perform dual averaging on dt
	double log_dt = mu - Stat * sqrt(i + 1) / (shr * (i + 1 + i0));
	dt = exp(log_dt);
	double eta = pow(i + 1, -kappa);
	log_dt_avg = eta * log_dt + (1 - eta) * log_dt_avg;
  }

  // set dt on its final dual averaged value
  dt = exp(log_dt_avg);

  delete[] en_i;
  delete[] en_f;
}

double Hamiltonian::run_HMC(const DiracOperator &D,
							const Action &A,
							const int &Nt,
							const double &dt,
							const int &iter,
							gsl_rng *engine,
							const string &integrator) const {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto *en_i = new double[4];
  auto *en_f = new double[4];

  // return statistic
  double Stat = 0;

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
	// if it's not the first interation set potential to
	// previous final value, otherwise compute it
	if (i) {
	  en_i[0] = en_f[0];
	  en_i[1] = en_f[1];
	} else {
	  en_i[0] = A.dirac2(D);
	  en_i[1] = A.dirac4(D);
	}

	// core part of HMC
	Stat += run_HMC_core(D, A, Nt, dt, engine, en_i, en_f, integrator);
  }

  delete[] en_i;
  delete[] en_f;

  return (Stat / iter);
}

double Hamiltonian::run_HMC(const DiracOperator &D,
							const Action &A,
							const int &Nt,
							const double &dt_min,
							const double &dt_max,
							const int &iter,
							gsl_rng *engine,
							const string &integrator) const {
  // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian
  auto *en_i = new double[4];
  auto *en_f = new double[4];

  // return statistic
  double Stat = 0;

  // iter repetitions of leapfrog
  for (int i = 0; i < iter; ++i) {
	// if it's not the first interation set potential to
	// previous final value, otherwise compute it
	if (i) {
	  en_i[0] = en_f[0];
	  en_i[1] = en_f[1];
	} else {
	  en_i[0] = A.dirac2(D);
	  en_i[1] = A.dirac4(D);
	}


	// core part of HMC
	Stat += run_HMC_core(D, A, Nt, dt_min, dt_max, engine, en_i, en_f, integrator);
  }

  delete[] en_i;
  delete[] en_f;

  return (Stat / iter);
}

double Hamiltonian::run_HMC_duav_core(const DiracOperator &D,
									  const Action &A,
									  const int &Nt,
									  const double &dt,
									  gsl_rng *engine,
									  double *en_i,
									  double *en_f,
									  const string &integrator) const {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sample_mom(D, engine);

  // store previous configuration
  auto *mat_bk = new cx_mat[D.nHL];
  auto *mat = D.get_mats();
  for (int j = 0; j < D.nHL; j++) {
	mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculate_K(D);
  en_i[3] = A.get_g2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (integrator == "leapfrog") {
	leapfrog(D, Nt, dt, A.get_g2());
  } else if (integrator == "omelyan") {
	omelyan(D, Nt, dt, A.get_g2());
  }

  // calculate final hamiltonian
  en_f[0] = A.dirac2(D);
  en_f[1] = A.dirac4(D);
  en_f[2] = calculate_K(D);
  en_f[3] = A.get_g2() * en_f[0] + en_f[1] + en_f[2];


  // metropolis test

  // sometimes leapfrog diverges and Hf becomes nan.
  // so first of all address this case
  if (std::isnan(en_f[3])) {
	e = 0;
	// restore old configuration
	for (int j = 0; j < D.nHL; ++j)
	  mat[j] = mat_bk[j];
	en_f[0] = en_i[0];
	en_f[1] = en_i[1];
	en_f[2] = en_i[2];
	en_f[3] = en_i[3];
  }
	// now do the standard metropolis test
  else if (en_f[3] > en_i[3]) {
	double r = gsl_rng_uniform(engine);
	e = exp(en_i[3] - en_f[3]);

	if (r > e) {
	  // restore old configuration
	  for (int j = 0; j < D.nHL; ++j)
		mat[j] = mat_bk[j];
	  en_f[0] = en_i[0];
	  en_f[1] = en_i[1];
	  en_f[2] = en_i[2];
	  en_f[3] = en_i[3];
	}
  }

  delete[] mat_bk;

  return e;
}

double Hamiltonian::run_HMC_core(const DiracOperator &D,
								 const Action &A,
								 const int &Nt,
								 const double &dt,
								 gsl_rng *engine,
								 double *en_i,
								 double *en_f,
								 const string &integrator) const  {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sample_mom(D, engine);

  // store previous configuration
  auto *mat_bk = new cx_mat[D.nHL];
  auto *mat = D.get_mats();
  for (int j = 0; j < D.nHL; j++) {
	mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculate_K(D);
  en_i[3] = A.get_g2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (integrator == "leapfrog") {
	leapfrog(D, Nt, dt, A.get_g2());
  } else if (integrator == "omelyan") {
	omelyan(D, Nt, dt, A.get_g2());
  }

  // calculate final hamiltonian
  en_f[0] = A.dirac2(D);
  en_f[1] = A.dirac4(D);
  en_f[2] = calculate_K(D);
  en_f[3] = A.get_g2() * en_f[0] + en_f[1] + en_f[2];


  // metropolis test
  if (en_f[3] > en_i[3]) {
	double r = gsl_rng_uniform(engine);
	e = exp(en_i[3] - en_f[3]);

	if (r > e) {
	  // restore old configuration
	  for (int j = 0; j < D.nHL; ++j)
		mat[j] = mat_bk[j];
	  en_f[0] = en_i[0];
	  en_f[1] = en_i[1];
	  en_f[2] = en_i[2];
	  en_f[3] = en_i[3];
	}
  }

  delete[] mat_bk;

  return e;
}

double Hamiltonian::run_HMC_core_debug(const DiracOperator &D,
									   const Action &A,
									   const int &Nt,
									   const double &dt,
									   gsl_rng *engine,
									   const string &integrator) const {
  // exp(-dH) (return value)
  double e;

  // resample momentum
  sample_mom(D, engine);

  // store previous configuration
  auto *mat_bk = new cx_mat[D.nHL];
  auto *mat = D.get_mats();
  for (int j = 0; j < D.nHL; j++) {
	mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  double Si = A.calculate_S(D);
  double Ki = calculate_K(D);
  double Hi = Si + Ki;

  // integration
  if (integrator == "leapfrog") {
	leapfrog(D, Nt, dt, A.get_g2());
  } else if (integrator == "omelyan") {
	omelyan(D, Nt, dt, A.get_g2());
  }

  // calculate final hamiltonian
  double Sf = A.calculate_S(D);
  double Kf = calculate_K(D);
  double Hf = Sf + Kf;

  e = exp(Hi - Hf);

  // metropolis test
  if (Hf > Hi) {
	double r = gsl_rng_uniform(engine);

	if (r > e) {
	  // restore old configuration
	  for (int j = 0; j < D.nHL; ++j)
		mat[j] = mat_bk[j];
	}
  }

  delete[] mat_bk;

  return e;
}

double Hamiltonian::run_HMC_core(const DiracOperator &D,
								 const Action &A,
								 const int &Nt,
								 const double &dt_min,
								 const double &dt_max,
								 gsl_rng *engine,
								 double *en_i,
								 double *en_f,
								 const string &integrator) const {
  // acceptance probability (return value)
  double e = 1;

  // resample momentum
  sample_mom(D, engine);

  // choose uniformly from [dt_min, dt_max)
  double dt = dt_min + (dt_max - dt_min) * gsl_rng_uniform(engine);

  // store previous configuration
  auto *mat_bk = new cx_mat[D.nHL];
  auto *mat = D.get_mats();
  for (int j = 0; j < D.nHL; j++) {
	mat_bk[j] = mat[j];
  }

  // calculate initial hamiltonian
  en_i[2] = calculate_K(D);
  en_i[3] = A.get_g2() * en_i[0] + en_i[1] + en_i[2];

  // integration
  if (integrator == "leapfrog") {
	leapfrog(D, Nt, dt, A.get_g2());
  } else if (integrator == "omelyan") {
	omelyan(D, Nt, dt, A.get_g2());
  }

  // calculate final hamiltonian
  en_f[0] = A.dirac2(D);
  en_f[1] = A.dirac4(D);
  en_f[2] = calculate_K(D);
  en_f[3] = A.get_g2() * en_f[0] + en_f[1] + en_f[2];


  // metropolis test
  if (en_f[3] > en_i[3]) {
	double r = gsl_rng_uniform(engine);
	e = exp(en_i[3] - en_f[3]);

	if (r > e) {
	  // restore old configuration
	  for (int j = 0; j < D.nHL; ++j)
		mat[j] = mat_bk[j];
	  en_f[0] = en_i[0];
	  en_f[1] = en_i[1];
	  en_f[2] = en_i[2];
	  en_f[3] = en_i[3];
	}
  }

  delete[] mat_bk;

  return e;
}

void Hamiltonian::setStepSize(double dt_) {this-> dt = dt_;}
void Hamiltonian::setIntegrator(Integrator integrator_) {this->integrator = integrator_;}
void Hamiltonian::setEngine(gsl_rng *engine_) { this->engine = engine_;}

