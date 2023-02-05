//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_HAMILTONIAN_HPP
#define RFL_HAMILTONIAN_HPP
#include "IAlgorithm.hpp"
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

enum Integrator {
  LEAPFROG,
  OMELYAN
};

class Hamiltonian : public IAlgorithm {
public:
  Hamiltonian(Integrator integrator, const gsl_rng* engine, double step_size);

  double updateDirac(const DiracOperator& dirac, const Action& action) const override;

  void setEngine(const gsl_rng* engine);
  const gsl_rng* getEngine() const { return this->m_engine; };

  void setIntegrator(Integrator integrator);
  Integrator getIntegrator() const { return this->m_integrator; };

  void setStepSize(double dt);
  double getStepSize() const { return this->m_dt; };

private:
  Integrator m_integrator = Integrator::LEAPFROG;
  const gsl_rng* m_engine;
  double m_dt;

  // This method seems to be the initialiser for the mom variables in DiracOperator
  void sampleMoments(const DiracOperator& dirac) const;
  double calculateK(const DiracOperator& dirac) const;
  double calculateH(const DiracOperator& dirac, const Action& action) const;

  double runHmc(const DiracOperator& dirac,
                const Action& action,
                const int& num_iterations,
                const int& iter) const;

  double runHmcDuavCore(const DiracOperator& dirac,
                        const Action& action,
                        const int& nt,
                        double* en_i,
                        double* en_f) const;

  double runHmcCore(const DiracOperator& dirac,
                    const Action& action,
                    const int& nt,
                    double* en_i,
                    double* en_f) const;

  double runHmcCoreDebug(const DiracOperator& dirac,
                         const Action& action,
                         const int& nt) const;

  // The methods below modify the step size "this->dt".
  void runHmcDuav(const DiracOperator& dirac,
                  const Action& action,
                  const int& nt,
                  const int& iter,
                  const double& target);

  double runHmc(const DiracOperator& dirac,
                const Action& action,
                const int& nt,
                const double& dt_min,
                const double& dt_max,
                const int& iter);

  double runHmcCore(const DiracOperator& dirac,
                    const Action& action,
                    const int& nt,
                    const double& dt_min,
                    const double& dt_max,
                    double* en_i,
                    double* en_f);

  // INTEGRATORS
  void leapfrog(const DiracOperator& dirac,
                const int& nt,
                double g_2) const;

  void omelyan(const DiracOperator& dirac,
               const int& nt,
               double g_2) const;
};

#endif//RFL_HAMILTONIAN_HPP