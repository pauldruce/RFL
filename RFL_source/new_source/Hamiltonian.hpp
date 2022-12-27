//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_HAMILTONIAN_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_HAMILTONIAN_HPP_
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "DiracOperator.hpp"
#include "Action.hpp"
//#include <string>

class Hamiltonian {
 public:
  Hamiltonian() = default;
 private:
  // This method seems to be the initialiser for the mom variables in DiracOperator
  void sample_mom(const DiracOperator &D,
				  gsl_rng *engine);
  double calculate_K(const DiracOperator &D);
  double calculate_H(const DiracOperator &D, const Action &A);

  void leapfrog(const DiracOperator &D,
				const int &Nt,
				const double &dt,
				const double g2);

  void omelyan(const DiracOperator &D,
			   const int &Nt,
			   const double &dt,
			   const double g2);

  void run_HMC_duav(const DiracOperator &D,
					const Action &A,
					const int &Nt,
					double &dt,
					const int &iter,
					gsl_rng *engine,
					const double &target,
					const std::string &integrator);

  double run_HMC(const DiracOperator &D,
				 const Action &A,
				 const int &Nt,
				 const double &dt,
				 const int &iter,
				 gsl_rng *engine,
				 const std::string &integrator);

  double run_HMC(const DiracOperator &D,
				 const Action &A,
				 const int &Nt,
				 const double &dt_min,
				 const double &dt_max,
				 const int &iter,
				 gsl_rng *engine,
				 const std::string &integrator);

  double run_HMC_duav_core(const DiracOperator &D,
						   const Action &A,
						   const int &Nt,
						   const double &dt,
						   gsl_rng *engine,
						   double *en_i,
						   double *en_f,
						   const std::string &integrator);

  double run_HMC_core(const DiracOperator &D,
					  const Action &A,
					  const int &Nt,
					  const double &dt,
					  gsl_rng *engine,
					  double *en_i,
					  double *en_f,
					  const std::string &integrator);

  double run_HMC_core_debug(const DiracOperator &D,
							const Action &A,
							const int &Nt,
							const double &dt,
							gsl_rng *engine,
							const std::string &integrator);

  double run_HMC_core(const DiracOperator &D,
					  const Action &A,
					  const int &Nt,
					  const double &dt_min,
					  const double &dt_max,
					  gsl_rng *engine,
					  double *en_i,
					  double *en_f,
					  const std::string &integrator);

};

#endif //RFL_RFL_SOURCE_NEW_SOURCE_HAMILTONIAN_HPP_
