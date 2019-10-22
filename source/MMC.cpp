#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "geometry.hpp"
#include "clifford.hpp"

using namespace std;
using namespace arma;

double Geom24::MMC_core(const double& scale, gsl_rng* engine, double* s_i, double* s_f)
{
    // acceptance probability
    double e;

    // metropolis
    int x = nHL*gsl_rng_uniform(engine);
    int I = dim*gsl_rng_uniform(engine);
    int J = dim*gsl_rng_uniform(engine);
    double re = 0;
    double im = 0;
    cx_double z;
    if(I != J)
    {
        re = scale*(-1.+2.*gsl_rng_uniform(engine));
        im = scale*(-1.+2.*gsl_rng_uniform(engine));
        z = cx_double(re, im);
    }
    else
    {
        re = scale*(-1.+2.*gsl_rng_uniform(engine));
        z = cx_double(re, 0);
    }

    double dS2 = delta2(x, I, J, z);
    double dS4 = delta4(x, I, J, z);
    double dS = g2*dS2+dS4;

    // metropolis test
    if(dS < 0)
    {
        // update matrix element
        if(I != J)
        {
            mat[x](I,J) += z;
            mat[x](J,I) += conj(z);
        }
        else
            mat[x](I,I) += 2.*z;

        // update action
        s_f[0] = s_i[0]+dS2;
        s_f[1] = s_i[1]+dS4;

        // move accepted
        e = 1;
    }
    else
    {
        e = exp(-dS);
        double p = gsl_rng_uniform(engine);

        if(e>p)
        {
            // update matrix element
            if(I != J)
            {
                mat[x](I,J) += z;
                mat[x](J,I) += conj(z);
            }
            else
                mat[x](I,I) += 2.*z;

            // update action
            s_f[0] = s_i[0]+dS2;
            s_f[1] = s_i[1]+dS4;
        }
        else
        {
            s_f[0] = s_i[0];
            s_f[1] = s_i[1];
        }
    }

    return e;
}
        


// MMC routine that performs dual averaging and outputs S2, S4, H, L
double Geom24::MMC(double& scale, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s, ostream& out_hl, const double& asymp, const double& shr/*=0.05*/, const double& kappa/*=0.75*/, const int& i0/*=10*/)
{
    // initial (_i) and final (_f) action2 and action4 
    double* s_i = new double [2];
    double* s_f = new double [2];

    // calculate length of a sweep in terms of dofs
    int Nsw = nHL*dim*dim;

    // dual averaging variables
    double Stat = 0;
    double mu = log(10*scale);
    double log_scale_avg = log(scale);

    // iter sweeps of metropolis
    for(int i=0; i<iter; ++i)
    {
        for(int j=0; j<Nsw; ++j)
        {
            // set action to previous final value,
            // unless it's the first iteration
            if(j)
            {
                s_i[0] = s_f[0];
                s_i[1] = s_f[1];
            }
            else
            {
                s_i[0] = dirac2();
                s_i[1] = dirac4();
            }

            
            Stat += asymp - MMC_core(scale, engine, s_i, s_f);
        
            // perform dual averaging
            double log_scale = mu - Stat*sqrt(i+1)/(shr*(i+1+i0));
            scale = exp(log_scale);
            double eta = pow(i+1, -kappa);
            log_scale_avg = eta*log_scale + (1-eta)*log_scale_avg;
        }

        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();

        // print every "gap" iteration
        if( !(i%gap) )
        {
            // action2 and action4
            out_s << s_f[0] << " " << s_f[1] << endl; 

            // mat
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
    
    // set scale on its final dual averaged value
    scale = exp(log_scale_avg);
    
    delete [] s_i;
    delete [] s_f;

    return (asymp - Stat/(iter*Nsw));
}

// MMC routine that performs dual averaging and outputs S2, S4
double Geom24::MMC(double& scale, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s, const double& asymp, const double& shr/*=0.05*/, const double& kappa/*=0.75*/, const int& i0/*=10*/)
{
    // initial (_i) and final (_f) action2 and action4 
    double* s_i = new double [2];
    double* s_f = new double [2];

    // calculate length of a sweep in terms of dofs
    int Nsw = nHL*dim*dim;

    // dual averaging variables
    double Stat = 0;
    double mu = log(10*scale);
    double log_scale_avg = log(scale);

    // iter sweeps of metropolis
    for(int i=0; i<iter; ++i)
    {
        for(int j=0; j<Nsw; ++j)
        {
            // set action to previous final value,
            // unless it's the first iteration
            if(j)
            {
                s_i[0] = s_f[0];
                s_i[1] = s_f[1];
            }
            else
            {
                s_i[0] = dirac2();
                s_i[1] = dirac4();
            }

            
            Stat += asymp - MMC_core(scale, engine, s_i, s_f);
        
            // perform dual averaging
            double log_scale = mu - Stat*sqrt(i+1)/(shr*(i+1+i0));
            scale = exp(log_scale);
            double eta = pow(i+1, -kappa);
            log_scale_avg = eta*log_scale + (1-eta)*log_scale_avg;
        }

        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();

        // print every "gap" iteration
        if( !(i%gap) )
        {
            // action2 and action4
            out_s << s_f[0] << " " << s_f[1] << endl; 
        }
    }
    
    // set scale on its final dual averaged value
    scale = exp(log_scale_avg);
    
    delete [] s_i;
    delete [] s_f;

    return (asymp - Stat/(iter*Nsw));
}

// MMC routine that performs dual averaging and doesn't output
double Geom24::MMC(double& scale, const int& iter, const int& adj, gsl_rng* engine, const double& asymp, const double& shr/*=0.05*/, const double& kappa/*=0.75*/, const int& i0/*=10*/)
{
    // initial (_i) and final (_f) action2 and action4 
    double* s_i = new double [2];
    double* s_f = new double [2];

    // calculate length of a sweep in terms of dofs
    int Nsw = nHL*dim*dim;

    // dual averaging variables
    double Stat = 0;
    double mu = log(10*scale);
    double log_scale_avg = log(scale);

    // iter sweeps of metropolis
    for(int i=0; i<iter; ++i)
    {
        for(int j=0; j<Nsw; ++j)
        {
            // set action to previous final value,
            // unless it's the first iteration
            if(j)
            {
                s_i[0] = s_f[0];
                s_i[1] = s_f[1];
            }
            else
            {
                s_i[0] = dirac2();
                s_i[1] = dirac4();
            }

            
            Stat += asymp - MMC_core(scale, engine, s_i, s_f);
        
            // perform dual averaging
            double log_scale = mu - Stat*sqrt(i+1)/(shr*(i+1+i0));
            scale = exp(log_scale);
            double eta = pow(i+1, -kappa);
            log_scale_avg = eta*log_scale + (1-eta)*log_scale_avg;
        }
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();
    }
    
    // set scale on its final dual averaged value
    scale = exp(log_scale_avg);
    
    delete [] s_i;
    delete [] s_f;

    return (asymp - Stat/(iter*Nsw));
}

// MMC routine that doesn't perform dual averaging and outputs S2, S4, H, L
double Geom24::MMC(const double& scale, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s, ostream& out_hl)
{
    // initial (_i) and final (_f) action2 and action4 
    double* s_i = new double [2];
    double* s_f = new double [2];

    // calculate length of a sweep in terms of dofs
    int Nsw = nHL*dim*dim;

    // return statistic
    double Stat = 0;

    // iter sweeps of metropolis
    for(int i=0; i<iter; ++i)
    {
        for(int j=0; j<Nsw; ++j)
        {
            // set action to previous final value,
            // unless it's the first iteration
            if(j)
            {
                s_i[0] = s_f[0];
                s_i[1] = s_f[1];
            }
            else
            {
                s_i[0] = dirac2();
                s_i[1] = dirac4();
            }

            
            Stat += MMC_core(scale, engine, s_i, s_f);
        }

        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();
        
        // print every "gap" iteration
        if( !(i%gap) )
        {
            // action2 and action4
            out_s << s_f[0] << " " << s_f[1] << endl; 

            // mat
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
    
    delete [] s_i;
    delete [] s_f;

    return (Stat/(iter*Nsw));
}

// MMC routine that doesn't perform dual averaging and outputs S2, S4
double Geom24::MMC(const double& scale, const int& iter, const int& gap, const int& adj, gsl_rng* engine, ostream& out_s)
{
    // initial (_i) and final (_f) action2 and action4 
    double* s_i = new double [2];
    double* s_f = new double [2];

    // calculate length of a sweep in terms of dofs
    int Nsw = nHL*dim*dim;

    // return statistic
    double Stat = 0;

    // iter sweeps of metropolis
    for(int i=0; i<iter; ++i)
    {
        for(int j=0; j<Nsw; ++j)
        {
            // set action to previous final value,
            // unless it's the first iteration
            if(j)
            {
                s_i[0] = s_f[0];
                s_i[1] = s_f[1];
            }
            else
            {
                s_i[0] = dirac2();
                s_i[1] = dirac4();
            }

            
            Stat += MMC_core(scale, engine, s_i, s_f);
        }

        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();
        
        // print every "gap" iteration
        if( !(i%gap) )
        {
            // action2 and action4
            out_s << s_f[0] << " " << s_f[1] << endl; 
        }
    }
    
    delete [] s_i;
    delete [] s_f;

    return (Stat/(iter*Nsw));
}

// MMC routine that doesn't perform dual averaging and doesn't output
double Geom24::MMC(const double& scale, const int& iter, const int& adj, gsl_rng* engine)
{
    // initial (_i) and final (_f) action2 and action4 
    double* s_i = new double [2];
    double* s_f = new double [2];

    // calculate length of a sweep in terms of dofs
    int Nsw = nHL*dim*dim;

    // return statistic
    double Stat = 0;

    // iter sweeps of metropolis
    for(int i=0; i<iter; ++i)
    {
        for(int j=0; j<Nsw; ++j)
        {
            // set action to previous final value,
            // unless it's the first iteration
            if(j)
            {
                s_i[0] = s_f[0];
                s_i[1] = s_f[1];
            }
            else
            {
                s_i[0] = dirac2();
                s_i[1] = dirac4();
            }

            
            Stat += MMC_core(scale, engine, s_i, s_f);
        }
        
        // adjust once every "adj" iterations
        if( !(i%adj) ) adjust();
    }
    
    delete [] s_i;
    delete [] s_f;

    return (Stat/(iter*Nsw));
}
