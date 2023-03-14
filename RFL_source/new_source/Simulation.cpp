//
// Created by Paul Druce on 07/12/2022.
//

#include "Simulation.hpp"

Simulation::Simulation(const DiracOperator& dirac, const Action& action, std::unique_ptr<IAlgorithm> &&monte_carlo_algorithm)
    : m_dirac(dirac), m_action(action), m_algorithm(std::move(monte_carlo_algorithm)) { };
