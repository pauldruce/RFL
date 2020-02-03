#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"

using namespace std;
using namespace arma;

#define N 10000
#define M 5

int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));

    // create geometry
    cout << "insert p, q, dim, g2" << endl;
    Geom24 G(cin);
    ifstream in("in.txt");
    G.read_mat(in);
    for(int i=0; i<G.get_nHL(); ++i)
        cout << G.get_mat(i) << endl;
    
    cout << G.der_dirac4(0, true) << endl;
            
    return 0;
}
