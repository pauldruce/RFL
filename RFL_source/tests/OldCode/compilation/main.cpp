#include "Geom24.hpp"
#include <iostream>

using namespace std;

int main() {
  int p, q, dim;
  double g2;
  cout << "Input p, q, dim, g2" << endl;
  cin >> p >> q >> dim >> g2;

  Geom24 G(p, q, dim, g2);
  cout << G << endl;

  return 0;
}
