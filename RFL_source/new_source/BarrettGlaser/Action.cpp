//
// Created by Paul Druce on 12/11/2022.
//
#include "Action.hpp"

using namespace std;
using namespace arma;

Action::Action(const double g_2) : m_g_2(g_2), m_g_4(1.0) {}

Action::Action(const double g_2, const double g_4) : m_g_2(g_2), m_g_4(g_4) {}

double Action::calculateSFromDirac(const IDiracOperator& dirac) const {
  const cx_mat dirac_mat = dirac.getDiracMatrix();
  const cx_mat dirac_squared = dirac_mat * dirac_mat;
  const double trace_dirac_squared = trace(dirac_squared).real();
  const double trace_dirac_4 = trace(dirac_squared * dirac_squared).real();
  return m_g_2 * trace_dirac_squared + m_g_4 * trace_dirac_4;
}

double Action::calculateS(const IDiracOperator& dirac) const {
  return m_g_2 * dirac.traceOfDiracSquared() + m_g_4 * dirac.traceOfDirac4();
}

void Action::setParams(const double g_2, const double g_4) {
  this->m_g_2 = g_2;
  this->m_g_4 = g_4;
}
void Action::setG4(const double value) { this->m_g_4 = value; }
void Action::setG2(const double value) { this->m_g_2 = value; }
