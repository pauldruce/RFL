//
// Created by Paul Druce on 17/11/2023.
//

#include "IDiracOperatorDerivatives.hpp"
#include "BarrettGlaser/BarrettGlaserAction.hpp"

using namespace arma;

cx_mat derDirac24(const IDiracOperator& dirac, const int& k, const bool& herm, const double g_2) {
  return BarrettGlaserAction::derDirac24(dirac, k, herm, g_2, 1.0);
}

arma::cx_mat derDirac2(const IDiracOperator& dirac, const int& k) {
  return BarrettGlaserAction::derDirac2(dirac, k);
}

cx_mat derDirac4(const IDiracOperator& dirac, const int& k, const bool& herm) {
  return BarrettGlaserAction::derDirac4(dirac, k, herm);
}
