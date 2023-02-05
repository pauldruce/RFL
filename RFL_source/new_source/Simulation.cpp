//
// Created by Paul Druce on 07/12/2022.
//

#include "Simulation.hpp"

Simulation::Simulation(const DiracOperator& D, const Action& A, IAlgorithm& M)
    : m_D(D), m_A(A), m_M(M) { };
