//
// Created by Paul Druce on 12/02/2023.
//

// For Mauro's implementation
#include <iostream>
#include <gsl/gsl_rng.h>
#include <ctime>
#include "Geom24.hpp"

// For the new implementation
#include "DiracOperator.hpp"
#include "Action.hpp"
#include "GslRng.hpp"
#include "Metropolis.hpp"
#include "Simulation.hpp"

// For benchmarking the two implementations
#include <chrono>

constexpr int P = 2;
constexpr int Q = 0;
constexpr int MATRIX_DIM = 25;
constexpr double G_2 = -3;
constexpr double SCALE = -3;
constexpr int NUM_STEPS = 10;
constexpr int NUM_ITERATIONS = 20;


void NewMethod(){
  auto rng = std::make_unique<GslRng>();

  DiracOperator dirac(P, Q,MATRIX_DIM);
  Action action(G_2, 1.0);
  auto metropolis = std::make_unique<Metropolis>(SCALE, NUM_STEPS, std::move(rng));
  auto simulation = Simulation(dirac,action,std::move(metropolis));

  for (int i = 0; i < NUM_ITERATIONS; i++) {
    simulation.run();
  }
}


void MauroMethod(){
  // Initialize the random number generator
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));
  // Create the Dirac operator
  Geom24 geom(P, Q, MATRIX_DIM, G_2);
  // Metropolis simulation
  for (int i = 0; i < NUM_ITERATIONS; ++i) {
    // Metropolis evolution for 100 steps
    geom.MMC(SCALE, NUM_STEPS, engine);
  }
}

int main(){
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  auto mauro_start = high_resolution_clock::now();
  MauroMethod();
  auto mauro_stop = high_resolution_clock::now();

  /* Getting number of milliseconds as an integer. */
  auto mauro_ms_int = duration_cast<milliseconds>(mauro_stop - mauro_start);

  std::cout << "Mauro's implementation:\n";
  std::cout << mauro_ms_int.count() << "ms\n";

  auto new_start = high_resolution_clock::now();
  NewMethod();
  auto new_stop = high_resolution_clock::now();

  auto new_ms_int = duration_cast<milliseconds>(new_stop-new_start);
  std::cout << "New implementation:\n";
  std::cout << new_ms_int.count() << "ms\n";

  return 0;
}