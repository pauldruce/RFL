//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_HAMILTONIAN_HPP
#define RFL_HAMILTONIAN_HPP
#include <memory>

#include "../IAlgorithm.hpp"
#include "../IRng.hpp"
#include "./Action.hpp"

enum Integrator {
  LEAPFROG,
  OMELYAN
};

/**
 * @class Hamiltonian
 *
 * @brief Implements the Hybrid Monte Carlo (HMC) algorithm for the Barrett-Glaser action.
 *
 * Uses molecular dynamics integration and Metropolis accept/reject steps to sample
 * Dirac operator configurations.
 */
class Hamiltonian final : public IAlgorithm {
public:
  /**
   * Default constructor deleted.
   */
  Hamiltonian() = delete;

  /**
   * Constructs a Hamiltonian simulation algorithm instance.
   *
   * @param action Unique pointer to a Barrett-Glaser Action object.
   * @param integrator Numerical integration scheme (LEAPFROG or OMELYAN).
   * @param step_size Integration step size dt.
   * @param rng Unique pointer to random number generator engine.
   */
  Hamiltonian(std::unique_ptr<Action>&& action, Integrator integrator, double step_size, std::unique_ptr<IRng>&& rng);

  /**
   * Runs HMC trajectory updates on the Dirac operator.
   *
   * @param dirac Dirac operator to update.
   * @return Mean acceptance rate across trajectories.
   */
  double updateDirac(const IDiracOperator& dirac) const override;

  /**
   * Sets the numerical integrator scheme.
   *
   * @param integrator Integrator scheme (LEAPFROG or OMELYAN).
   */
  void setIntegrator(Integrator integrator);

  /**
   * Returns the current numerical integrator scheme.
   *
   * @return Current Integrator enum value.
   */
  Integrator getIntegrator() const { return this->m_integrator; };

  /**
   * Sets the integration step size dt.
   *
   * @param dt Integration step size.
   */
  void setStepSize(double dt);

  /**
   * Returns the current integration step size dt.
   *
   * @return Integration step size.
   */
  double getStepSize() const { return this->m_dt; };

private:
  std::unique_ptr<Action> m_action;
  Integrator m_integrator = LEAPFROG;
  double m_dt;
  std::unique_ptr<IRng> m_rng;

  // Samples Gaussian conjugate momenta for the Dirac operator.
  void sampleMoments(const IDiracOperator& dirac) const;
  static double calculateK(const IDiracOperator& dirac);
  double calculateH(const IDiracOperator& dirac) const;

  double run(const IDiracOperator& dirac,
             const int& num_iterations,
             const int& iter) const;

  double runDualAveragingCore(const IDiracOperator& dirac,
                              const int& nt,
                              std::vector<double>& en_i,
                              std::vector<double>& en_f) const;

  double runCore(const IDiracOperator& dirac,
                 const int& nt,
                 std::vector<double>& en_i,
                 std::vector<double>& en_f) const;

  double runCoreDebug(const IDiracOperator& dirac,
                      const int& nt) const;

  // The methods below modify the step size m_dt.
  void runDualAverage(const IDiracOperator& dirac,
                      const int& nt,
                      const int& iter,
                      const double& target);

  double run(const IDiracOperator& dirac,
             const int& nt,
             const double& dt_min,
             const double& dt_max,
             const int& iter);

  double runCore(const IDiracOperator& dirac,
                 const int& nt,
                 const double& dt_min,
                 const double& dt_max,
                 std::vector<double>& en_i,
                 std::vector<double>& en_f);

  // Numerical integrators.
  void leapfrog(const IDiracOperator& dirac,
                const int& nt,
                double g_2) const;

  void omelyan(const IDiracOperator& dirac,
               const int& nt,
               double g_2) const;
};

#endif//RFL_HAMILTONIAN_HPP
