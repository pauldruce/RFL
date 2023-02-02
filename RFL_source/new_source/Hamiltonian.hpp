//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_RFL_SOURCE_NEW_SOURCE_HAMILTONIAN_HPP_
#define RFL_RFL_SOURCE_NEW_SOURCE_HAMILTONIAN_HPP_
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "IAlgorithm.hpp"

enum Integrator {
  leapfrog,
  omelyan
} ;

class Hamiltonian : public IAlgorithm {
 public:
  Hamiltonian() = default;
  double updateDirac(const DiracOperator &D, const Action &A) const override {
	// TODO: Change occurrences of string "leapfrog" etc, to enums.
	std::string str_integrator = "leapfrog";
	if(integrator != Integrator::leapfrog){
	  str_integrator = "omelyan";
	}
	const double acceptance_val_per_iter = this->run_HMC(
		D,A,
		10,
		dt,
		1000,
		engine,
		str_integrator
	);
	return acceptance_val_per_iter;
  };

  void setEngine(gsl_rng *engine);
  void setIntegrator(Integrator integrator);
  void setStepSize(double step_size);

 private:
  gsl_rng *engine;
  Integrator integrator = Integrator::leapfrog;
  double dt;

  // This method seems to be the initialiser for the mom variables in DiracOperator
  void sample_mom(const DiracOperator &D,
				  gsl_rng *engine) const;
  double calculate_K(const DiracOperator &D) const;
  double calculate_H(const DiracOperator &D, const Action &A) const;

  void leapfrog(const DiracOperator &D,
				const int &Nt,
				const double &dt,
				const double g2) const;

  void omelyan(const DiracOperator &D,
			   const int &Nt,
			   const double &dt,
			   const double g2) const;

  void run_HMC_duav(const DiracOperator &D,
					const Action &A,
					const int &Nt,
					double &dt,
					const int &iter,
					gsl_rng *engine,
					const double &target,
					const std::string &integrator) const;

  double run_HMC(const DiracOperator &D,
				 const Action &A,
				 const int &Nt,
				 const double &dt,
				 const int &iter,
				 gsl_rng *engine,
				 const std::string &integrator) const;

  double run_HMC(const DiracOperator &D,
				 const Action &A,
				 const int &Nt,
				 const double &dt_min,
				 const double &dt_max,
				 const int &iter,
				 gsl_rng *engine,
				 const std::string &integrator) const;

  double run_HMC_duav_core(const DiracOperator &D,
						   const Action &A,
						   const int &Nt,
						   const double &dt,
						   gsl_rng *engine,
						   double *en_i,
						   double *en_f,
						   const std::string &integrator) const;

  double run_HMC_core(const DiracOperator &D,
					  const Action &A,
					  const int &Nt,
					  const double &dt,
					  gsl_rng *engine,
					  double *en_i,
					  double *en_f,
					  const std::string &integrator) const;

  double run_HMC_core_debug(const DiracOperator &D,
							const Action &A,
							const int &Nt,
							const double &dt,
							gsl_rng *engine,
							const std::string &integrator) const;

  double run_HMC_core(const DiracOperator &D,
					  const Action &A,
					  const int &Nt,
					  const double &dt_min,
					  const double &dt_max,
					  gsl_rng *engine,
					  double *en_i,
					  double *en_f,
					  const std::string &integrator) const ;

};

#endif //RFL_RFL_SOURCE_NEW_SOURCE_HAMILTONIAN_HPP_
