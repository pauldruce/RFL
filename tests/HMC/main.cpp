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

#define N 100000
#define M 5

int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));

    // create geometry
    cout << "insert p, q, dim, g2" << endl;
    Geom24 G(cin);
    G.shuffle(engine);
        
    // thermalize first
    double tgt = 0.8;
    double dt = 0.000001;
    G.HMC_duav(10, dt, 10000, engine, tgt);
    cout << "dual averaging complete" << endl;
    cout << "dual averaged dt: " << dt << endl;
    double ar = G.HMC(10, dt, 10000, engine);
    cout << "thermalization complete" << endl;
    cout << "acceptance rate: " << ar << endl;
    
    vec O(N);
    for(int j=0; j<N; ++j)
    {
        double o = G.HMC_core_debug(10, dt, engine);
        
        O(j) = o;
    }
    cout << "avg exp(-dH): " << mean(O) << endl;

    return 0;
}
