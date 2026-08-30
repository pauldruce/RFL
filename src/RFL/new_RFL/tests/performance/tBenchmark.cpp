//
// Created by Paul Druce on 01/04/2023.
//

#include <gtest/gtest.h>

// Legacy reference implementation
#include "Geom24.hpp"
#include <ctime>
#include <gsl/gsl_rng.h>
#include <iostream>

// New RFL implementation
#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"
#include "Simulation.hpp"

// Benchmark timing utilities
#include <chrono>

constexpr int P = 2;
constexpr int Q = 0;
constexpr int MATRIX_DIM = 25;
constexpr double G_2 = -3;
constexpr double SCALE = -3;
constexpr int NUM_STEPS = 10;
constexpr int NUM_ITERATIONS = 5;
constexpr double MAX_FRACTION_DIFF = 1.05;

void NewMethod() {
  auto rng = std::make_unique<GslRng>();
  auto dirac = std::make_unique<DiracOperator>(P, Q, MATRIX_DIM);
  auto action = std::make_unique<Action>(G_2, 1.0);
  auto metropolis = std::make_unique<Metropolis>(std::move(action), SCALE, NUM_STEPS, std::move(rng));

  const auto simulation = Simulation(
      std::move(dirac),
      std::move(metropolis));

  for (int i = 0; i < NUM_ITERATIONS; i++) {
    simulation.run();
  }
}

void MauroMethod() {
  // Initialise the random number generator
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));
  // Create the Dirac operator
  Geom24 geom(P, Q, MATRIX_DIM, G_2);
  // Run Metropolis sampling
  for (int i = 0; i < NUM_ITERATIONS; ++i) {
    // Run Metropolis sweep steps
    geom.MMC(SCALE, NUM_STEPS, engine);
  }

  gsl_rng_free(engine);
}

TEST(BenchmarkTests, NewImplementationIsNotSignificantlySlower) {
  using std::chrono::duration;
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::milliseconds;

  const auto mauro_start = high_resolution_clock::now();
  MauroMethod();
  const auto mauro_stop = high_resolution_clock::now();

  // Calculate elapsed duration in milliseconds.
  const auto mauro_ms_int = duration_cast<milliseconds>(mauro_stop - mauro_start);

  std::cout << "Mauro's implementation:\n";
  std::cout << mauro_ms_int.count() << "ms\n";

  const auto new_start = high_resolution_clock::now();
  NewMethod();
  const auto new_stop = high_resolution_clock::now();

  const auto new_ms_int = duration_cast<milliseconds>(new_stop - new_start);
  std::cout << "New implementation:\n";
  std::cout << new_ms_int.count() << "ms\n";

  const auto winner = new_ms_int.count() < mauro_ms_int.count() ? "the new implementation" : "Mauro's implementation";
  std::cout << "The winner is: " << winner << std::endl;

  ASSERT_TRUE(new_ms_int.count() < mauro_ms_int.count() * MAX_FRACTION_DIFF) << "New implementation is slower by more than 5 percent.";
}