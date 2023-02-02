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
  Hamiltonian(Integrator integrator, gsl_rng* engine, double step_size);

  double updateDirac(const DiracOperator &D, const Action &A) const override;;

  void setEngine(gsl_rng *engine);
  gsl_rng* getEngine() const { return this->engine; };

  void setIntegrator(Integrator integrator);
  Integrator getIntegrator() const { return this->integrator; };

  void setStepSize(double step_size);
  double getStepSize() const { return this->dt; };

 private:
  Integrator integrator = Integrator::leapfrog;
  gsl_rng *engine;
  double dt;

  // This method seems to be the initialiser for the mom variables in DiracOperator
  void sample_mom(const DiracOperator &D) const;
  double calculate_K(const DiracOperator &D) const;
  double calculate_H(const DiracOperator &D, const Action &A) const;

  double run_HMC(const DiracOperator &D,
				 const Action &A,
				 const int &Nt,
				 const int &iter) const;

  double run_HMC_duav_core(const DiracOperator &D,
						   const Action &A,
						   const int &Nt,
						   double *en_i,
						   double *en_f) const;

  double run_HMC_core(const DiracOperator &D,
					  const Action &A,
					  const int &Nt,
					  double *en_i,
					  double *en_f) const;

  double run_HMC_core_debug(const DiracOperator &D,
							const Action &A,
							const int &Nt) const;

  // The methods below modify the step size "this->dt".
  void run_HMC_duav(const DiracOperator &D,
					const Action &A,
					const int &Nt,
					const int &iter,
					const double &target);

  double run_HMC(const DiracOperator &D,
				 const Action &A,
				 const int &Nt,
				 const double &dt_min,
				 const double &dt_max,
				 const int &iter);

  double run_HMC_core(const DiracOperator &D,
					  const Action &A,
					  const int &Nt,
					  const double &dt_min,
					  const double &dt_max,
					  double *en_i,
					  double *en_f);


  // INTEGRATORS
  void leapfrog(const DiracOperator &D,
				const int &Nt,
				double g2) const;

  void omelyan(const DiracOperator &D,
			   const int &Nt,
			   double g2) const;
};

#endif //RFL_RFL_SOURCE_NEW_SOURCE_HAMILTONIAN_HPP_
