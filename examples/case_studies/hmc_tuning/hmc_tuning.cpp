//
// Created by Paul Druce on 09/12/2022.
//
#include "Geom24.hpp"
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  // Initialise the random number generator.
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));

  Geom24 G(2, 0, 10, -2.7);

  // Determine output directory from environment or default to current directory.
  const char* out_env = std::getenv("RFL_OUTPUT_DIR");
  std::filesystem::path out_dir = (out_env && *out_env) ? std::filesystem::path(out_env) : std::filesystem::current_path();
  if (!std::filesystem::exists(out_dir)) {
    std::filesystem::create_directories(out_dir);
  }

  // Open output files.
  ofstream out_S(out_dir / "example_S.txt");
  ofstream out_HL(out_dir / "example_HL.txt");

  // Tune step scale with dual-averaging.
  double tgt = 0.8; // Target acceptance rate.
  double dt = 0.001;// Initial guess for dt.
  G.HMC_duav(10, dt, 10000, engine, tgt, "leapfrog");
  cout << "dual-averaging complete" << endl;
  cout << "dual-averaged dt: " << dt << endl;
  // Thermalisation.
  double acc_rate = G.HMC(10, dt, 10000, engine, "leapfrog");
  cout << "thermalisation complete" << endl;
  cout << "acceptance rate: " << acc_rate << endl;
  // Hamiltonian Monte Carlo simulation.
  for (int i = 1; i < 1000; ++i) {
    G.HMC(10, dt, 1000, engine, "leapfrog");
    G.print_S(out_S);
    G.print_HL(out_HL);
  }
  out_S.close();
  out_HL.close();
}
