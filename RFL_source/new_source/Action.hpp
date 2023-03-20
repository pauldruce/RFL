//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_ACTION_HPP
#define RFL_ACTION_HPP
#include "IAction.hpp"
#include "DiracOperator.hpp"
#include <armadillo>
#include <string>

/**
 * Action implements the \f$S(D) = g_2Tr(D^2) + g_4 Tr(D^4) \f$ action as described in the papers
 * of John Barrett and L Glaser.
 */
class Action : public IAction {
public:
  // CONSTRUCTORS AND DESTRUCTORS

  Action(double g_2, double g_4);
  explicit Action(double g_2);
  Action() : m_g_2(0.0), m_g_4(0.0) { };
  ~Action() = default;

  // METHODS
  /**
   * setG2 updates the value of the quadratic coupling variables - which is typically called
   * \f$g_2\f$ in the literature.
   * @param value The value of g2 you want to set
   */
  void setG2(double value);
  void setG4(double value);
  void setParams(double g_2, double g_4);
  double getG2() const { return m_g_2; }
  double getG4() const { return m_g_4; }

  double calculateS(const DiracOperator& dirac) const override;
  double calculateSFromDirac(const DiracOperator& dirac) const;

  void printS(const DiracOperator& dirac, std::ostream& out) const {
    out << dirac.traceOfDiracSquared() << " " << dirac.traceOfDirac4() << std::endl;
  }

private:
  double m_g_2, m_g_4;

  double computeA4(const DiracOperator& dirac,
                   const int& i_1,
                   const int& i_2,
                   const int& i_3,
                   const int& i_4) const;
  double computeA2(const DiracOperator& dirac, const int& i_1, const int& i_2) const;
  double computeA(const DiracOperator& dirac, const int& i) const;
};
#endif// RFL_ACTION_HPP
