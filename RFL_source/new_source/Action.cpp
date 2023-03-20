//
// Created by Paul Druce on 12/11/2022.
//
#include "Action.hpp"

using namespace std;
using namespace arma;

Action::Action(double g_2) : m_g_2(g_2), m_g_4(1.0) { }

Action::Action(double g_2, double g_4) : m_g_2(g_2), m_g_4(g_4) { }

double Action::calculateSFromDirac(const DiracOperator& dirac) const {
  cx_mat dirac_mat = dirac.getDiracMatrix();
  cx_mat dirac_squared = dirac_mat * dirac_mat;
  double trace_dirac_squared = trace(dirac_squared).real();
  double trace_dirac_4 = trace(dirac_squared * dirac_squared).real();
  return m_g_2 * trace_dirac_squared + m_g_4 * trace_dirac_4;
}

double Action::calculateS(const DiracOperator& dirac) const {
  return m_g_2 * dirac.traceOfDiracSquared() + m_g_4 * dirac.traceOfDirac4();
}

void Action::setParams(double g_2, double g_4) {
  this->m_g_2 = g_2;
  this->m_g_4 = g_4;
}
void Action::setG4(double value) { this->m_g_4 = value; }
void Action::setG2(double value) { this->m_g_2 = value; }
