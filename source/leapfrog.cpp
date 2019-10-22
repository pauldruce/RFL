#include "geometry.hpp"

using namespace std;
using namespace arma;

//  Standard leapfrog:
//  K/2
//  for (Nt-1) times:
//    S
//    K
//  S
//  K/2
void Geom24::leapfrog(const int& Nt, const double& dt)
{
    for(int i=0; i<nHL; ++i)
        mat[i] += (dt/2.)*mom[i].st();

    for(int j=0; j<Nt-1; ++j)
    {
        for(int i=0; i<nHL; ++i)
        {
            mom[i] += -dt*der_dirac24(i, true);
            //mom[i] += -dt*der_dirac2(i);
            mat[i] += dt*mom[i].st();
        }
    }
    
    for(int i=0; i<nHL; ++i)
    {
        mom[i] += -dt*der_dirac24(i, true);
        //mom[i] += -dt*der_dirac2(i);
        mat[i] += (dt/2.)*mom[i].st();
    }
}


// Split leapfrog to expolit cheap potential.
// In the following, S4=U0 and g2*S2=U1.
// This is the splitting based on Neal 5.1:
// U1/2
// for (Nt-1) times:
//   U0/2M
//   for (M-1) times:
//     K/M
//     U0/M
//   K/M
//   U0/2M
//   U1
// U0/2M
// for (M-1) times:
//   K/M
//   U0/M
// K/M
// U0/2M
// U1/2
void Geom24::leapfrog(const int& Nt, const double& dt, const int& M)
{

    for(int i=0; i<nHL; ++i)
        mom[i] += (-dt/2)*der_dirac2(i);

    for(int j=0; j<Nt-1; ++j)
    {
        for(int i=0; i<nHL; ++i)
            mom[i] += (-dt*g2/(2*M))*der_dirac4(i, true);

        for(int k=0; k<M-1; ++k)
        {
            for(int i=0; i<nHL; ++i)
            {
                mat[i] += (dt/M)*mom[i].st();
                mom[i] += (-dt*g2/M)*der_dirac4(i, true);
            }
        }
            
        for(int i=0; i<nHL; ++i)
        {
            mat[i] += (dt/M)*mom[i].st();
            mom[i] += (-dt*g2/(2*M))*der_dirac4(i, true);
            mom[i] += -dt*der_dirac2(i);
        }
    }

    for(int i=0; i<nHL; ++i)
        mom[i] += (-dt*g2/(2*M))*der_dirac4(i, true);

    for(int k=0; k<M-1; ++k)
    {
        for(int i1=0; i1<nHL; ++i1)
        {
            mat[i1] += (dt/M)*mom[i1].st();
            mom[i1] += (-dt*g2/M)*der_dirac4(i1, true);
        }
    }
        
    for(int i=0; i<nHL; ++i)
    {
        mat[i] += (dt/M)*mom[i].st();
        mom[i] += (-dt*g2/(2*M))*der_dirac4(i, true);
        mom[i] += (-dt/2)*der_dirac2(i);
    }
}
