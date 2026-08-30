#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>

#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "Clifford.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"

using namespace std;
using namespace arma;

int main() {
  cout << "========================================" << endl;
  cout << " RFL C++ API Example" << endl;
  cout << "========================================" << endl
       << endl;

  // Configure safety limits
  int default_max = Clifford::getMaxMode();
  cout << "Default Max Clifford Mode (p+q): " << default_max << endl;

  int new_max = 20;
  cout << "Overriding Max Clifford Mode to: " << new_max << " (allows larger matrices)\n"
       << endl;
  Clifford::setMaxMode(new_max);

  // Initialize a Dirac Operator
  int p = 1;
  int q = 3;
  int dim = 10;
  cout << "Initializing Dirac Operator with p=" << p << ", q=" << q << ", dim=" << dim << "..." << endl;

  auto dirac = make_unique<DiracOperator>(p, q, dim);

  cout << "-> Matrix Dimension: " << dirac->getMatrixDimension() << "\n"
       << endl;

  // Setup the Metropolis Algorithm
  double g_2 = -1.0;
  double g_4 = 1.0;
  double scale = 1.0;
  int num_steps = 100;
  unsigned long seed = 42;

  cout << "Running Metropolis Algorithm (g_2=" << g_2 << ", g_4=" << g_4 << ", steps=" << num_steps << ")..." << endl;

  auto action = make_unique<Action>(g_2, g_4);
  auto rng = make_unique<GslRng>(seed);
  Metropolis metropolis(move(action), scale, num_steps, move(rng));

  double acceptance_rate = metropolis.updateDirac(*dirac);
  cout << "-> Acceptance Rate: " << fixed << setprecision(2) << (acceptance_rate * 100.0) << "%\n"
       << endl;

  // Analyze Eigenvalues
  cout << "Computing eigenvalues..." << endl;
  vec eigenvals_real = real(dirac->getEigenvalues());

  cout << "-> Found " << eigenvals_real.n_elem << " eigenvalues." << endl;
  cout << "-> Min eigenvalue: " << fixed << setprecision(4) << eigenvals_real.min() << endl;
  cout << "-> Max eigenvalue: " << fixed << setprecision(4) << eigenvals_real.max() << endl;
  cout << "-> Mean eigenvalue: " << fixed << setprecision(4) << mean(eigenvals_real) << endl;
  cout << "\nSimulation complete!" << endl;

  return 0;
}
