//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_HAMILTONIAN_HPP
#define RFL_HAMILTONIAN_HPP
#include "IAlgorithm.hpp"
#include "IRng.hpp"
#include "Action.hpp"
#include <cmath>
#include <memory>

enum Integrator {
  LEAPFROG,
  OMELYAN
};

class Hamiltonian : public IAlgorithm {
public:
  Hamiltonian() = delete;
  Hamiltonian(Integrator integrator, std::unique_ptr<IRng> &&rng, double step_size);

  /**
   * The updateDirac method for Hamiltonian is not implemented for all classes that satisfy the
   * IAction interface. This is only currently implemented for the class Action.
   *
   * TODO: Generalise the hamiltonian method to accept other options.
   *
   * @param dirac a reference to the DiracOperator to be updated.
   * @param action a reference to the IAction-derived class.
   * @return
   */
  double updateDirac(const DiracOperator& dirac, const IAction& action) const override {
      return this->updateDirac(
          dirac,
          dynamic_cast<const Action&>(action)
      );
  }

  double updateDirac(const DiracOperator& dirac, const Action& action) const;

  void setIntegrator(Integrator integrator);
  Integrator getIntegrator() const { return this->m_integrator; };

  void setStepSize(double dt);
  double getStepSize() const { return this->m_dt; };

private:
  Integrator m_integrator = Integrator::LEAPFROG;
  std::unique_ptr<IRng> m_rng;
  double m_dt;

  // This method seems to be the initialiser for the mom variables in DiracOperator
  void sampleMoments(const DiracOperator& dirac) const;
  double calculateK(const DiracOperator& dirac) const;
  double calculateH(const DiracOperator& dirac, const Action& action) const;

  double run(const DiracOperator& dirac,
             const Action& action,
             const int& num_iterations,
             const int& iter) const;

  double runDualAveragingCore(const DiracOperator& dirac,
                              const Action& action,
                              const int& nt,
                              double* en_i,
                              double* en_f) const;

  double runCore(const DiracOperator& dirac,
                 const Action& action,
                 const int& nt,
                 double* en_i,
                 double* en_f) const;

  double runCoreDebug(const DiracOperator& dirac,
                      const Action& action,
                      const int& nt) const;

  // The methods below modify the step size "this->dt".
  void runDualAverage(const DiracOperator& dirac,
                      const Action& action,
                      const int& nt,
                      const int& iter,
                      const double& target);

  double run(const DiracOperator& dirac,
             const Action& action,
             const int& nt,
             const double& dt_min,
             const double& dt_max,
             const int& iter);

  double runCore(const DiracOperator& dirac,
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
