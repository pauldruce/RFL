//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_ACTION_HPP
#define RFL_ACTION_HPP
#include "IAction.hpp"
#include "IDiracOperator.hpp"

/**
 * Action implements the \f$ S(D) = g_2Tr(D^2) + g_4 Tr(D^4) \f$ action as described in the papers
 * of John Barrett and L Glaser.
 */
class Action : public IAction {
public:
  // CONSTRUCTORS AND DESTRUCTORS

  /**
   * The preferred constructor for the Barrett-Glaser action class.
   *
   * @param g_2 The value of the quadratic coupling constant, \f$g_2\f$ you want to set.
   * @param g_4 The value of the quartic coupling constant, \f$g_4\f$ you want to set.
   */
  Action(double g_2, double g_4);

  /**
   * A constructor for the Barrett-Glaser action class that assumes the quartic coupling constant is equal to 1.
   *
   * @param g_2 The value of the quadratic coupling constant, \f$g_2\f$ you want to set.
   */
  explicit Action(double g_2);

  /**
   * A constructor for the Barrett-Glaser action that sets the coupling constants to equal zero.
   */
  Action() : m_g_2(0.0), m_g_4(0.0){};

  ~Action() override = default;

  // METHODS
  /**
   * setG2 updates the value of the quadratic coupling variables - which is typically called
   * \f$g_2\f$ in the literature.
   * @param value The value of \f$g_2\f$ you want to set
   */
  void setG2(double value);

  /**
   * setG4 updates the value of the quartic coupling variables - which is typically called
   * \f$g_4\f$ in the literature.
   * @param value The value of \f$g_4\f$ you want to set
   */
  void setG4(double value);

  /**
   * setParams updates the values of the two coupling constants in the Barrett-Glaser action,
   * @param g_2 The value of the quadratic coupling constant, \f$g_2\f$ you want to set.
   * @param g_4 The value of the quartic coupling constant, \f$g_4\f$ you want to set.
   */
  void setParams(double g_2, double g_4);

  /**
   * getG2 returns the current value of the quadratic coupling constant of the action
   */
  double getG2() const { return m_g_2; }

  /**
   * getG4 returns the current value of the quartic coupling constant of the action.
   */
  double getG4() const { return m_g_4; }

  /**
   * This calculates the Barrett-Glaser action for the supplied Dirac operator. This makes use of an
   * an optimised algorithm for the B-G action.
   * @param dirac A reference to the Dirac operator you want to calculate the Barrett-Glaser action of.
   */
  double calculateS(const IDiracOperator& dirac) const override;

  /**
   * This calculates the Barrett-Glaser action for the supplied Dirac operator by computing the full
   * matrix representation of the Dirac operator and directly calculating the action. This is not an
   * optimum calculation. Please make use of calculateS(...) above for faster results.
   * @param dirac A reference to the Dirac operator you want to calculate the Barrett-Glaser action of.
   */
  double calculateSFromDirac(const IDiracOperator& dirac) const;

private:
  double m_g_2, m_g_4;
};
#endif// RFL_ACTION_HPP
