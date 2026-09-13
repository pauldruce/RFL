//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_ACTION_HPP
#define RFL_ACTION_HPP
#include "../IAction.hpp"
#include "../IDiracOperator.hpp"

/**
 * @class Action
 *
 * @brief Implements the Barrett-Glaser action \f$ S(D) = g_2\text{Tr}(D^2) + g_4 \text{Tr}(D^4) \f$.
 */
class Action : public IAction {
public:
  /**
   * Constructs a Barrett-Glaser action with quadratic and quartic coupling constants.
   *
   * @param g_2 Quadratic coupling constant.
   * @param g_4 Quartic coupling constant.
   */
  Action(double g_2, double g_4);

  /**
   * Constructs a Barrett-Glaser action with quartic coupling constant equal to 1.0.
   *
   * @param g_2 Quadratic coupling constant.
   */
  explicit Action(double g_2);

  /**
   * Constructs a Barrett-Glaser action with coupling constants set to 0.0.
   */
  Action() : m_g_2(0.0), m_g_4(0.0){};

  ~Action() override = default;

  /**
   * Sets the quadratic coupling constant \f$g_2\f$.
   *
   * @param value Quadratic coupling constant.
   */
  void setG2(double value);

  /**
   * Sets the quartic coupling constant \f$g_4\f$.
   *
   * @param value Quartic coupling constant.
   */
  void setG4(double value);

  /**
   * Sets the quadratic and quartic coupling constants.
   *
   * @param g_2 Quadratic coupling constant.
   * @param g_4 Quartic coupling constant.
   */
  void setParams(double g_2, double g_4);

  /**
   * Returns the current quadratic coupling constant \f$g_2\f$.
   *
   * @return Quadratic coupling constant.
   */
  double getG2() const { return m_g_2; }

  /**
   * Returns the current quartic coupling constant \f$g_4\f$.
   *
   * @return Quartic coupling constant.
   */
  double getG4() const { return m_g_4; }

  /**
   * Calculates the Barrett-Glaser action using component matrix traces.
   *
   * @param dirac Dirac operator reference.
   * @return Action value.
   */
  double calculateS(const IDiracOperator& dirac) const override;

  /**
   * Calculates the Barrett-Glaser action from the assembled Dirac operator matrix.
   *
   * Use calculateS() for faster performance.
   *
   * @param dirac Dirac operator reference.
   * @return Action value.
   */
  double calculateSFromDirac(const IDiracOperator& dirac) const;

private:
  double m_g_2, m_g_4;
};
#endif// RFL_ACTION_HPP
