#include <gsl/gsl_rng.h>
#include "geometry.hpp"

using namespace std;
using namespace arma;


// HMC routine that doesn't performs dual averaging and outputs S2, S4, H, L
double Geom24::HMC_rand(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s, ostream& out_hl)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // return statistic
    double Stat = 0;
    
    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        // if it's not the first interation set potential to
        // previous final value, otherwise compute it
        if(i)
        {
            en_i[0] = en_f[0];
            en_i[1] = en_f[1];
        }
        else
        {
            en_i[0] = dirac2();
            en_i[1] = dirac4();
        }

        
        // core part of HMC
        Stat += HMC_rand_nosplit_core(Nt_min, Nt_max, dt_min, dt_max, engine, en_i, en_f);
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();

        // print once every "gap" iterations
        if( !(i%gap) )
        {
            // print S2 and S4
            out_s << en_f[0] << " " << en_f[1] << endl;
            
            // print mat
            for(int j=0; j<nHL; ++j)
            {
                for(int k=0; k<dim; ++k)
                {
                    for(int l=0; l<dim; ++l)
                        out_hl << mat[j](k,l).real() << " " << mat[j](k,l).imag() << " ";
                }
                out_hl << endl;
            }
        }
    }

    delete [] en_i;
    delete [] en_f;

    return (Stat/iter);
}

// HMC routine that doesn't performs dual averaging and outputs S2, S4
double Geom24::HMC_rand(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // return statistic
    double Stat = 0;
    
    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        // if it's not the first interation set potential to
        // previous final value, otherwise compute it
        if(i)
        {
            en_i[0] = en_f[0];
            en_i[1] = en_f[1];
        }
        else
        {
            en_i[0] = dirac2();
            en_i[1] = dirac4();
        }

        
        // core part of HMC
        Stat += HMC_rand_nosplit_core(Nt_min, Nt_max, dt_min, dt_max, engine, en_i, en_f);
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();

        // print once every "gap" iterations
        if( !(i%gap) )
        {
            // print S2 and S4
            out_s << en_f[0] << " " << en_f[1] << endl;
        }
    }

    delete [] en_i;
    delete [] en_f;

    return (Stat/iter);
}

// HMC routine that doesn't performs dual averaging and doesn't output
double Geom24::HMC_rand(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, const int& iter, const int& adj, gsl_rng* engine)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // return statistic
    double Stat = 0;
    
    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        // if it's not the first interation set potential to
        // previous final value, otherwise compute it
        if(i)
        {
            en_i[0] = en_f[0];
            en_i[1] = en_f[1];
        }
        else
        {
            en_i[0] = dirac2();
            en_i[1] = dirac4();
        }

        
        // core part of HMC
        Stat += HMC_rand_nosplit_core(Nt_min, Nt_max, dt_min, dt_max, engine, en_i, en_f);
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();
    }

    delete [] en_i;
    delete [] en_f;

    return (Stat/iter);
}
