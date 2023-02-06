//
// Created by Paul Druce on 07/12/2022.
//

#include "Simulation.hpp"

Simulation::Simulation(const DiracOperator& dirac, const Action& action, IAlgorithm& monte_carlo_algorithm)
    : m_dirac(dirac), m_action(action), m_algorithm(monte_carlo_algorithm) { };
