#include <gsl/gsl_rng.h>
#include "geometry.hpp"

using namespace std;
using namespace arma;


// HMC routine that performs dual averaging and outputs S2, S4, H, L
double Geom24::HMC_fix(const int& Nt, double& dt, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s, ostream& out_hl, const double& asymp, const double& shr/*=0.05*/, const double& kappa/*=0.75*/, const int& i0/*=10*/)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // dual averaging variables for dt
    double Stat = 0;
    double mu = log(10*dt);
    double log_dt_avg = log(dt);
    
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
        Stat += asymp - HMC_fix_nosplit_core(Nt, dt, engine, en_i, en_f);
        
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

        // perform dual averaging on dt
        double log_dt = mu - Stat*sqrt(i+1)/(shr*(i+1+i0));
        dt = exp(log_dt);
        double eta = pow(i+1, -kappa);
        log_dt_avg = eta*log_dt + (1-eta)*log_dt_avg;
    }

    // set dt on its final dual averaged value
    dt = exp(log_dt_avg);


    delete [] en_i;
    delete [] en_f;

    return (asymp - Stat/iter);
}

// HMC routine that performs dual averaging and outputs S2, S4
double Geom24::HMC_fix(const int& Nt, double& dt, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s, const double& asymp, const double& shr/*=0.05*/, const double& kappa/*=0.75*/, const int& i0/*=10*/)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // dual averaging variables for dt
    double Stat = 0;
    double mu = log(10*dt);
    double log_dt_avg = log(dt);
    
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
        Stat += asymp - HMC_fix_nosplit_core(Nt, dt, engine, en_i, en_f);
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();

        // print once every "gap" iterations
        if( !(i%gap) )
        {
            // print S2 and S4
            out_s << en_f[0] << " " << en_f[1] << endl;
        }

        // perform dual averaging on dt
        double log_dt = mu - Stat*sqrt(i+1)/(shr*(i+1+i0));
        dt = exp(log_dt);
        double eta = pow(i+1, -kappa);
        log_dt_avg = eta*log_dt + (1-eta)*log_dt_avg;
    }

    // set dt on its final dual averaged value
    dt = exp(log_dt_avg);


    delete [] en_i;
    delete [] en_f;

    return (asymp - Stat/iter);
}

// HMC routine that performs dual averaging and doesn't output
double Geom24::HMC_fix(const int& Nt, double& dt, const int& iter, const int& adj, gsl_rng* engine, const double& asymp, const double& shr/*=0.05*/, const double& kappa/*=0.75*/, const int& i0/*=10*/)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // dual averaging variables for dt
    double Stat = 0;
    double mu = log(10*dt);
    double log_dt_avg = log(dt);
    
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
        Stat += asymp - HMC_fix_nosplit_core(Nt, dt, engine, en_i, en_f);
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();

        // perform dual averaging on dt
        double log_dt = mu - Stat*sqrt(i+1)/(shr*(i+1+i0));
        dt = exp(log_dt);
        double eta = pow(i+1, -kappa);
        log_dt_avg = eta*log_dt + (1-eta)*log_dt_avg;
    }

    // set dt on its final dual averaged value
    dt = exp(log_dt_avg);


    delete [] en_i;
    delete [] en_f;

    return (asymp - Stat/iter);
}

// HMC routine that doesn't performs dual averaging and outputs S2, S4, H, L
double Geom24::HMC_fix(const int& Nt, const double& dt, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s, ostream& out_hl)
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
        Stat += HMC_fix_nosplit_core(Nt, dt, engine, en_i, en_f);
        
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
double Geom24::HMC_fix(const int& Nt, const double& dt, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s)
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
        Stat += HMC_fix_nosplit_core(Nt, dt, engine, en_i, en_f);
        
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
double Geom24::HMC_fix(const int& Nt, const double& dt, const int& iter, const int& adj, gsl_rng* engine)
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
        Stat += HMC_fix_nosplit_core(Nt, dt, engine, en_i, en_f);
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();
    }

    delete [] en_i;
    delete [] en_f;

    return (Stat/iter);
}
