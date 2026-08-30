//
// Created by Paul Druce on 17/11/2023.
//

#ifndef IDIRACOPERATORDERIVATIVES_HPP
#define IDIRACOPERATORDERIVATIVES_HPP
#include "IDiracOperator.hpp"

/**
 * Calculates the derivative of the Barrett-Glaser action with respect to matrix element k.
 *
 * @param dirac Dirac operator reference.
 * @param k Index of the matrix.
 * @param herm True if the matrix is Hermitian, false otherwise.
 * @param g_2 Coupling constant for Tr(D^2).
 * @return Derivative matrix.
 */
arma::cx_mat derDirac24(const IDiracOperator& dirac, const int& k, const bool& herm, double g_2);

/**
 * Calculates the derivative of Tr(D^2) with respect to matrix element k.
 *
 * @param dirac Dirac operator reference.
 * @param k Index of the matrix.
 * @return Derivative matrix.
 */
arma::cx_mat derDirac2(const IDiracOperator& dirac, const int& k);

/**
 * Calculates the derivative of Tr(D^4) with respect to matrix element k.
 *
 * @param dirac Dirac operator reference.
 * @param k Index of the matrix.
 * @param herm True if the matrix is Hermitian, false otherwise.
 * @return Derivative matrix.
 */
arma::cx_mat derDirac4(const IDiracOperator& dirac, const int& k, const bool& herm);

#endif//IDIRACOPERATORDERIVATIVES_HPP
