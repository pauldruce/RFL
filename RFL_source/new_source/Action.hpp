//
// Created by Paul Druce on 12/11/2022.
//

#ifndef RFL_ACTION_HPP
#define RFL_ACTION_HPP
#include "DiracOperator.hpp"
#include <armadillo>
#include <gsl/gsl_rng.h>
#include <string>

class Action {
public:
  // CONSTRUCTORS AND DESTRUCTORS
  Action(double g_2, double g_4);
  explicit Action(double g_2);
  Action() : m_g2(0.0), m_g4(0.0) { };
  ~Action() = default;

  // METHODS
  void setG2(double value);
  void setG4(double value);
  void setParams(double g_2, double g_4);
  double getG2() const { return m_g2; }
  double getG4() const { return m_g4; }

  double calculateS(const DiracOperator& dirac) const;
  double calculateSFromDirac(const DiracOperator& dirac) const;
  double dirac2(const DiracOperator& dirac) const;
  double dirac4(const DiracOperator& dirac) const;

  void printS(const DiracOperator& dirac, std::ostream& out) const {
    out << dirac2(dirac) << " " << dirac4(dirac) << std::endl;
  }

private:
  double m_g2, m_g4;

  double computeA4(const DiracOperator& dirac,
                   const int& i_1,
                   const int& i_2,
                   const int& i_3,
                   const int& i_4) const;
  double computeA2(const DiracOperator& dirac, const int& i_1, const int& i_2) const;
  double computeA(const DiracOperator& dirac, const int& i) const;
};
#endif// RFL_ACTION_HPP
