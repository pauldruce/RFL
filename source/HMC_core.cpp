#include <gsl/gsl_rng.h>
#include "geometry.hpp"

using namespace std;
using namespace arma;


double Geom24::HMC_fix_core(const int& Nt, const double& dt, gsl_rng* engine, double* en_i, double* en_f)
{
    // acceptance probability (return value)
    double e;
    
    // resample momentum
    sample_mom(engine);

    // store previous configuration
    cx_mat* mat_bk = new cx_mat [nHL];
    for(int j=0; j<nHL; j++)
        mat_bk[j] = mat[j];

    // calculate initial hamiltonian
    en_i[2] = calculate_K();
    en_i[3] = g2*en_i[0]+en_i[1]+en_i[2];

    // leapfrog
    leapfrog(Nt, dt);

    // calculate final hamiltonian
    en_f[0] = dirac2();
    en_f[1] = dirac4();
    en_f[2] = calculate_K();
    en_f[3] = g2*en_f[0]+en_f[1]+en_f[2];


    // metropolis test
    
    // sometimes leapfrog diverges and Hf becomes nan.
    // so first of all address this case
    if(std::isnan(en_f[3]))
    {
        e = 0;
        // restore old configuration
        for(int j=0; j<nHL; ++j)
        {
            mat[j] = mat_bk[j];
            en_f[0] = en_i[0];
            en_f[1] = en_i[1];
            en_f[2] = en_i[2];
            en_f[3] = en_i[3];
        }
    }
    // now do the standard metropolis test
    else if(en_f[3] > en_i[3])
    {
        double r = gsl_rng_uniform(engine);
        e = exp(en_i[3]-en_f[3]);

        if(r > e)
        {
            // restore old configuration
            for(int j=0; j<nHL; ++j)
            {
                mat[j] = mat_bk[j];
                en_f[0] = en_i[0];
                en_f[1] = en_i[1];
                en_f[2] = en_i[2];
                en_f[3] = en_i[3];
            }
        }
    }
    else
        e = 1;

    delete [] mat_bk;

    return e;
}


double Geom24::HMC_rand_core(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, gsl_rng* engine, double* en_i, double* en_f)
{
    // acceptance probability (return value)
    double e;
    
    // resample momentum
    sample_mom(engine);

    // choose uniformly from [dt_min, dt_max)
    double dt = dt_min + (dt_max-dt_min)*gsl_rng_uniform(engine);

    // choose uniformly from [Nt_min, Nt_max)
    double Nt = Nt_min + (Nt_max-Nt_min)*gsl_rng_uniform(engine);

    // store previous configuration
    cx_mat* mat_bk = new cx_mat [nHL];
    for(int j=0; j<nHL; j++)
        mat_bk[j] = mat[j];

    // calculate initial hamiltonian
    en_i[2] = calculate_K();
    en_i[3] = g2*en_i[0]+en_i[1]+en_i[2];

    // leapfrog
    leapfrog(Nt, dt);

    // calculate final hamiltonian
    en_f[0] = dirac2();
    en_f[1] = dirac4();
    en_f[2] = calculate_K();
    en_f[3] = g2*en_f[0]+en_f[1]+en_f[2];


    // metropolis test
    
    // sometimes leapfrog diverges and Hf becomes nan.
    // so first of all address this case
    if(std::isnan(en_f[3]))
    {
        e = 0;
        // restore old configuration
        for(int j=0; j<nHL; ++j)
        {
            mat[j] = mat_bk[j];
            en_f[0] = en_i[0];
            en_f[1] = en_i[1];
            en_f[2] = en_i[2];
            en_f[3] = en_i[3];
        }
    }
    // now do the standard metropolis test
    else if(en_f[3] > en_i[3])
    {
        double r = gsl_rng_uniform(engine);
        e = exp(en_i[3]-en_f[3]);

        if(r > e)
        {
            // restore old configuration
            for(int j=0; j<nHL; ++j)
            {
                mat[j] = mat_bk[j];
                en_f[0] = en_i[0];
                en_f[1] = en_i[1];
                en_f[2] = en_i[2];
                en_f[3] = en_i[3];
            }
        }
    }
    else
        e = 1;

    delete [] mat_bk;

    return e;
}

