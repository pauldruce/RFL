//
// Barrett-Glaser action functional, local analytic variations, and gradients.
//

#ifndef RFL_BARRETT_GLASER_ACTION_HPP
#define RFL_BARRETT_GLASER_ACTION_HPP

#include "../IAction.hpp"
#include "../IDiracOperator.hpp"
#include <armadillo>

/**
 * @class BarrettGlaserAction
 *
 * @brief Evaluates the Barrett-Glaser action: `S(D) = g_2 Tr(D^2) + g_4 Tr(D^4)`.
 *
 * Computes the spectral action functional and its local variations for
 * Monte Carlo simulations of finite noncommutative geometries.
 *
 * @note Invariants:
 *       - Reality: S(D) is strictly real for Hermitian D.
 *       - Gauge Invariance: S(U D U†) == S(D) for unitary U.
 * @note Complexity: O(V N^3) for trace calculation, avoiding O((d_gamma N)^3) assembly.
 *
 * @see Literature: [Barrett2016, Eq. 1.1] (arXiv:1510.01377)
 * @see Theory: docs/theory/barrett_glaser.md#2-the-spectral-action-functional
 */
class BarrettGlaserAction : public IAction {
public:
  /**
   * @brief Constructs a Barrett-Glaser action with quadratic and quartic coupling constants.
   *
   * @param[in] g_2 Quadratic coupling constant.
   * @param[in] g_4 Quartic coupling constant.
   */
  BarrettGlaserAction(double g_2, double g_4);

  /**
   * @brief Constructs a Barrett-Glaser action with quartic coupling constant equal to 1.0.
   *
   * @param[in] g_2 Quadratic coupling constant.
   */
  explicit BarrettGlaserAction(double g_2);

  /**
   * @brief Constructs a Barrett-Glaser action with coupling constants initialised to 0.0.
   */
  BarrettGlaserAction() : m_g_2(0.0), m_g_4(0.0) {}

  ~BarrettGlaserAction() override = default;

  /**
   * @brief Sets the quadratic coupling constant g_2.
   *
   * @param[in] value Quadratic coupling constant.
   */
  void setG2(double value);

  /**
   * @brief Sets the quartic coupling constant g_4.
   *
   * @param[in] value Quartic coupling constant.
   */
  void setG4(double value);

  /**
   * @brief Sets the quadratic and quartic coupling constants.
   *
   * @param[in] g_2 Quadratic coupling constant.
   * @param[in] g_4 Quartic coupling constant.
   */
  void setParams(double g_2, double g_4);

  /**
   * @brief Returns the current quadratic coupling constant g_2.
   *
   * @return Quadratic coupling constant.
   */
  double getG2() const { return m_g_2; }

  /**
   * @brief Returns the current quartic coupling constant g_4.
   *
   * @return Quartic coupling constant.
   */
  double getG4() const { return m_g_4; }

  /**
   * @brief Calculates the Barrett-Glaser action using component matrix traces.
   *
   * Evaluates `S(D) = g_2 Tr(D^2) + g_4 Tr(D^4)` using Clifford trace factorisations
   * across the tensor product space without assembling the full Dirac operator.
   *
   * @param[in] dirac Dirac operator reference.
   * @return Calculated action value.
   *
   * @note Complexity: O(V N^3) vs O((d_gamma N)^3) for full Dirac assembly.
   *
   * @see Literature: [Barrett2016, Eq. 1.1] (arXiv:1510.01377)
   * @see Theory: docs/theory/barrett_glaser.md#2-the-spectral-action-functional
   */
  double calculateS(const IDiracOperator& dirac) const override;

  /**
   * @brief Calculates the Barrett-Glaser action from the assembled Dirac operator matrix.
   *
   * Assembles the full `(d_gamma N) x (d_gamma N)` Dirac operator matrix and evaluates
   * `S(D) = g_2 Re(Tr(D^2)) + g_4 Re(Tr(D^4))`. Used as a brute-force baseline for verification.
   *
   * @param[in] dirac Dirac operator reference.
   * @return Calculated action value.
   *
   * @note Complexity: O((d_gamma N)^3) matrix powers.
   *
   * @see Theory: docs/theory/barrett_glaser.md#2-the-spectral-action-functional
   */
  double calculateSFromDirac(const IDiracOperator& dirac) const;

  /**
   * @brief Calculates the analytic change in Tr(D^2) under an elementary Hermitian matrix variation.
   *
   * Evaluates `Delta Tr(D^2) = Tr((D + delta D)^2) - Tr(D^2)` under perturbation
   * `delta M = z E_ij + conj(z) E_ji`.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] x Matrix index (0 <= x < V).
   * @param[in] row_index Matrix row index i (0 <= row_index < N).
   * @param[in] column_index Matrix column index j (0 <= column_index < N).
   * @param[in] z Complex perturbation value.
   * @return Exact change in Tr(D^2).
   *
   * @note Complexity: O(1) for off-diagonal moves (i != j), O(N) for diagonal moves (i == j).
   *
   * @see Literature: [Barrett2016, Section 3, Eq. 3.3] (arXiv:1510.01377)
   * @see Theory: docs/theory/barrett_glaser.md#32-quadratic-variation-delta-s_2
   */
  static double delta2(
      const IDiracOperator& dirac,
      int x,
      int row_index,
      int column_index,
      const arma::cx_double& z);

  /**
   * @brief Calculates the analytic change in Tr(D^4) under an elementary Hermitian matrix variation.
   *
   * Evaluates `Delta Tr(D^4) = Tr((D + delta D)^4) - Tr(D^4)` under perturbation
   * `delta D = gamma^x (x) (z E_ij + conj(z) E_ji)`.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] x Matrix index (0 <= x < V).
   * @param[in] row_index Matrix row index i (0 <= row_index < N).
   * @param[in] column_index Matrix column index j (0 <= column_index < N).
   * @param[in] z Complex perturbation value.
   * @return Exact change in Tr(D^4).
   *
   * @note Complexity: O(V^3 N^3) on-the-fly matrix products per move vs O((d_gamma N)^3) Dirac recalculation. Can be optimised to O(V^3 N^2) with sweep-level matrix product caching.
   *
   * @see Literature: [Barrett2016, Section 3, Eq. 3.4-3.14] (arXiv:1510.01377)
   * @see Theory: docs/theory/barrett_glaser.md#33-quartic-variation-delta-s_4
   */
  static double delta4(
      const IDiracOperator& dirac,
      int x,
      int row_index,
      int column_index,
      const arma::cx_double& z);

  /**
   * @brief Calculates the total change in action given pre-computed quadratic and quartic variations.
   *
   * Linearly combines `Delta S = g_2 * delta_2 + g_4 * delta_4`.
   *
   * @param[in] delta_2 Pre-computed change in Tr(D^2).
   * @param[in] delta_4 Pre-computed change in Tr(D^4).
   * @return Total action change: g_2 * delta_2 + g_4 * delta_4.
   *
   * @see Theory: docs/theory/barrett_glaser.md#34-total-variation--acceptance
   */
  double delta24(double delta_2, double delta_4) const {
    return m_g_2 * delta_2 + m_g_4 * delta_4;
  }

  /**
   * @brief Calculates the total change in Barrett-Glaser action with explicit couplings.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] x Matrix index (0 <= x < V).
   * @param[in] row_index Matrix row index i (0 <= row_index < N).
   * @param[in] column_index Matrix column index j (0 <= column_index < N).
   * @param[in] z Complex perturbation value.
   * @param[in] g_2 Quadratic coupling constant.
   * @param[in] g_4 Quartic coupling constant.
   * @return Total action change: g_2 * delta_2 + g_4 * delta_4.
   *
   * @see Theory: docs/theory/barrett_glaser.md#34-total-variation--acceptance
   */
  static double delta24(
      const IDiracOperator& dirac,
      int x,
      int row_index,
      int column_index,
      const arma::cx_double& z,
      double g_2,
      double g_4);

  /**
   * @brief Calculates the total change in Barrett-Glaser action using instance couplings.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] x Matrix index (0 <= x < V).
   * @param[in] row_index Matrix row index i (0 <= row_index < N).
   * @param[in] column_index Matrix column index j (0 <= column_index < N).
   * @param[in] z Complex perturbation value.
   * @return Total action change: g_2 * delta_2 + g_4 * delta_4.
   *
   * @see Theory: docs/theory/barrett_glaser.md#34-total-variation--acceptance
   */
  double delta24(
      const IDiracOperator& dirac,
      int x,
      int row_index,
      int column_index,
      const arma::cx_double& z) const;

  /**
   * @brief Calculates the matrix derivative of Tr(D^2) with respect to matrix M_k.
   *
   * Evaluates `nabla_{M_k} Tr(D^2) = 4 d_gamma (N M_k + eps_k Tr(M_k) 1_N)`.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] k Matrix index (0 <= k < V).
   * @return Complex derivative matrix.
   *
   * @see Theory: docs/theory/barrett_glaser.md#41-quadratic-gradient
   */
  static arma::cx_mat derDirac2(const IDiracOperator& dirac, int k);

  /**
   * @brief Calculates the matrix derivative of Tr(D^4) with respect to matrix M_k.
   *
   * Evaluates `nabla_{M_k} Tr(D^4) = B_4(k) + B_2(k) + B(k)`. When `herm` is true,
   * projects onto the Hermitian tangent space: `nabla_{Herm} = 2 (nabla + nabla†)`.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] k Matrix index (0 <= k < V).
   * @param[in] herm True to project onto Hermitian tangent space, false otherwise.
   * @return Complex derivative matrix.
   *
   * @see Theory: docs/theory/barrett_glaser.md#42-quartic-gradient
   */
  static arma::cx_mat derDirac4(const IDiracOperator& dirac, int k, bool herm);

  /**
   * @brief Calculates the total action gradient with explicit couplings.
   *
   * Evaluates `nabla_{M_k} S(D) = g_2 nabla_{M_k} Tr(D^2) + g_4 nabla_{M_k} Tr(D^4)`.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] k Matrix index (0 <= k < V).
   * @param[in] herm True to project onto Hermitian tangent space, false otherwise.
   * @param[in] g_2 Quadratic coupling constant.
   * @param[in] g_4 Quartic coupling constant.
   * @return Complex gradient matrix.
   *
   * @see Theory: docs/theory/barrett_glaser.md#4-analytical-gradients--hmc-forces
   */
  static arma::cx_mat derDirac24(
      const IDiracOperator& dirac,
      int k,
      bool herm,
      double g_2,
      double g_4);

  /**
   * @brief Calculates the total action gradient using instance couplings.
   *
   * @param[in] dirac Dirac operator reference.
   * @param[in] k Matrix index (0 <= k < V).
   * @param[in] herm True to project onto Hermitian tangent space, false otherwise.
   * @return Complex gradient matrix.
   *
   * @see Theory: docs/theory/barrett_glaser.md#4-analytical-gradients--hmc-forces
   */
  arma::cx_mat derDirac24(const IDiracOperator& dirac, int k, bool herm) const;

private:
  double m_g_2;
  double m_g_4;
};

#endif// RFL_BARRETT_GLASER_ACTION_HPP
