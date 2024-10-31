//
// Created by Paul Druce on 09/12/2022.
//
#include "Geom24.hpp"
#include <ctime>
#include <gsl/gsl_rng.h>
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  // Initialize the random number generator
  gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(engine, time(nullptr));

  Geom24 G(2, 0, 10, -2.7);

  // Open output files
  ofstream out_S("example_S.txt");
  ofstream out_HL("example_HL.txt");

  // Tuning with dual averaging
  double tgt = 0.8; // Target acceptance rate
  double dt = 0.001;// Initial guess for dt
  G.HMC_duav(10, dt, 10000, engine, tgt, "leapfrog");
  cout << "dual averaging complete" << endl;
  cout << "dual averaged dt: " << dt << endl;
  // Thermalization
  double acc_rate = G.HMC(10, dt, 10000, engine, "leapfrog");
  cout << "thermalization complete" << endl;
  cout << "acceptance rate: " << acc_rate << endl;
  // Hamiltonian Monte Carlo simulation
  for (int i = 1; i < 1000; ++i) {
    G.HMC(10, dt, 1000, engine, "leapfrog");
    G.print_S(out_S);
    G.print_HL(out_HL);
  }
  out_S.close();
  out_HL.close();
}