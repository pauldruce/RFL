#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <iostream>
#include <string>

// Struct that packs the simulation
// parameters all together (hmc and mmc).
// 
// Explanation of integration parameters:
//
// L    =   # of integration steps in leapfrog
// 
// dL   =   randomize L in interval [L-dL, L+dL]
//
// dt   =   integration stepsize in leapfrog
//
// ddt  =   randomize dt in interval [dt-ddt, dt+ddt]
//
// AR   =   target acceptance rate
// 
// dAR  =   randomize dt in interval [dt_min, dt_max]
//          such that dt_min gives AR+dAR and
//          dt_max gives AR-dAR
//
// M    =   split hamiltonian in M iterations of
//          S2 for each iteration of S4
//
// Explanation of other parameters:
//
// scale        =   metropolis scale factor
//
// iter_therm   =   # of thermalization iterations
// 
// iter_simul   =   # of data-collecting iterations
//
// gap          =   # of iterations to be skipped between
//                  two measurements
//
// adj          =   # of iterations to be skipped between
//                  two applications of hermitization + tracelessization
//
// Explanation of MC mode:
//
// fix_nosplit  =   don't randomize dt or L, don't split hamiltonian
//
// fix_split    =   don't randomize dt or L, split hamiltonian
//
// rand_nosplit =   randomize dt and L, don't split hamiltonian
//
// rand_split   =   randomize dt and L, split hamiltonian
//
// mmc          =   metropolis
//
struct Simul_params
{
    // Geometric parameters
    int p;
    int q;
    int dim;

    // Monte Carlo parameters
    double AR;
    double dAR=0;

    // Hamiltonian parameters
    int L;
    int dL=0;
    double dt;
    double ddt=0;

    // Metropolis parameters
    double scale;

    // Other parameters
    int iter_therm;
    int iter_simul;
    int gap=1000;
    int adj=1000;

    // MC mode
    std::string mode;

    // Control string
    std::string control;
};

// Function to read simulation parameters from stream
bool read_init_stream(std::istream&, struct Simul_params&);

// Checks that the necessary parameters are there
bool params_validity(const struct Simul_params&);

#endif

