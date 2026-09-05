//
// Created by Paul Druce on 07/12/2022.
//

#include "Simulation.hpp"

// Initialise the simulation with a Dirac operator and a Monte Carlo algorithm.
Simulation::Simulation(std::unique_ptr<DiracOperator>&& dirac, std::unique_ptr<IAlgorithm>&& monte_carlo_algorithm)
    : m_dirac(std::move(dirac)), m_algorithm(std::move(monte_carlo_algorithm)){};
