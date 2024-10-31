//
// Created by Paul Druce on 17/11/2023.
//

#ifndef IDIRACOPERATORDERIVATIVES_HPP
#define IDIRACOPERATORDERIVATIVES_HPP
#include "IDiracOperator.hpp"

arma::cx_mat derDirac24(const IDiracOperator& dirac, const int& k, const bool& herm, double g_2);

arma::cx_mat derDirac2(const IDiracOperator& dirac, const int& k);

arma::cx_mat derDirac4(const IDiracOperator& dirac, const int& k, const bool& herm);

#endif//IDIRACOPERATORDERIVATIVES_HPP
