//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_HAMILTONIAN_HPP
#define RFL_HAMILTONIAN_HPP
#include <memory>

#include "Action.hpp"
#include "IAlgorithm.hpp"
#include "IRng.hpp"

enum Integrator {
  LEAPFROG,
  OMELYAN
};

/**
 * @class Hamiltonian
 * @brief Represents a Hamiltonian system for performing simulation.
 *
 * The Hamiltonian class encapsulates the necessary components for performing simulation
 * using Hamiltonian dynamics. It combines an action, an integrator, a step size, and a random
 * number generator into a coherent simulation algorithm.
 */
class Hamiltonian final : public IAlgorithm {
public:
  /**
   * The default constructor for this class has been disabled, please use another constructor.
   */
  Hamiltonian() = delete;

  /**
   * @class Hamiltonian
   * @brief Represents a Hamiltonian system for performing simulation.
   *
   * The Hamiltonian class encapsulates the necessary components for performing simulation
   * using Hamiltonian dynamics. It combines an action, an integrator, a step size, and a random
   * number generator into a coherent simulation algorithm.
   *
   * To construct this class, you need to pass in:
   * @param action a unique_ptr to the Barrett-Glaser action, Action
   * @param integrator an enum to indicate which Hamiltonian algorithm method you want to use.
   * @param step_size is the number of steps to take for each call of the method updateDirac.
   * @param rng a unique_ptr to a class that is derived from the abstract class IRng.
   */
  Hamiltonian(std::unique_ptr<Action>&& action, Integrator integrator, double step_size, std::unique_ptr<IRng>&& rng);

  double updateDirac(const IDiracOperator& dirac) const override;

  /**
 * @brief Sets the integrator for the system.
 *
 * This function sets the integrator for the system. The integrator is responsible
 * for computing the derivatives of the system variables and updating their values
 * over time.
 *
 * @param integrator The integrator object to be set.
 *
 * @sa Integrator
 */
  void setIntegrator(Integrator integrator);

  /**
   * @brief Get the integrator value.
   *
   * This function returns the current value of the integrator.
   *
   * @return The current value of the integrator.
   */
  Integrator getIntegrator() const { return this->m_integrator; };

  /**
   * @brief Sets the step size for the simulation.
   *
   * This function sets the step size, represented by the parameter `dt`, for the simulation.
   * The step size determines how different the next DiracOperator in the chain is from the
   * starting point.
   *
   * @param dt The step size for the simulation.
   */
  void setStepSize(double dt);

  /**
   * @brief Retrieves the step size value.
   *
   * This function returns the step size value associated with the current
   * object.
   *
   * @return The step size value.
   */
  double getStepSize() const { return this->m_dt; };

private:
  std::unique_ptr<Action> m_action;
  Integrator m_integrator = LEAPFROG;
  double m_dt;
  std::unique_ptr<IRng> m_rng;

  // This method seems to be the initialiser for the mom variables in DiracOperator.
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

  // The methods below modify the step size "this->dt".
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

  // INTEGRATORS
  void leapfrog(const IDiracOperator& dirac,
                const int& nt,
                double g_2) const;

  void omelyan(const IDiracOperator& dirac,
               const int& nt,
               double g_2) const;
};

#endif//RFL_HAMILTONIAN_HPP
