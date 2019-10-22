#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "geometry.hpp"
#include "clifford.hpp"

using namespace std;
using namespace arma;

// Constructors

Geom24::Geom24(int p_, int q_, int dim_, double g2_)
    : p(p_), q(q_), dim(dim_), g2(g2_)
{
    // initialize derived parameters
    derived_parameters();

    // allocate and initialize H and L matrices to identity
    mat = new arma::cx_mat [nHL];
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    for(int i=0; i<nHL; i++)
    {
        if(i<nH)
            eps[i] = 1;
        else
            eps[i] = -1;

        mat[i].eye(dim, dim); 
        mom[i].eye(dim, dim); 
    }


    //clog << "Geometry initialized with the following parameters:" << endl;
    //clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;
}

Geom24::Geom24(istream& in)
{
    read_parameters(in);
    
    // initialize derived parameters
    derived_parameters();

    // allocate and initialize H and L matrices to identity
    mat = new arma::cx_mat [nHL];
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    for(int i=0; i<nHL; i++)
    {
        if(i<nH)
            eps[i] = 1;
        else
            eps[i] = -1;

        mat[i].eye(dim, dim); 
        mom[i].eye(dim, dim); 
    }
    
    //clog << "Geometry initialized with the following parameters:" << endl;
    //clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;

}

// Copy constructor
Geom24::Geom24(const Geom24& G)
{
    // copy parameters
    dim = G.get_dim();
    p = G.get_p();
    q = G.get_q();
    g2 = G.get_g2();
    nH = G.get_nH();
    nL = G.get_nL();
    nHL = G.get_nHL();
    dim_omega = G.get_dim_omega();


    // allocate and copy matrices
    mat = new arma::cx_mat [nHL];
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    omega = new arma::cx_mat [nHL];
    for(int i=0; i<nHL; i++)
    {
        mat[i] = G.get_mat(i);
        mom[i] = G.get_mom(i);
        eps[i] = G.get_eps(i);
        omega[i] = G.get_omega(i);
    }
    
    int nHL4 = pow(nHL, 4);
    omega_table_4 = new cx_double [nHL4];
    for(int i=0; i<nHL4; i++)
        omega_table_4[i] = G.get_omega_table_4(i);

    
    //clog << "Geometry initialized with the following parameters:" << endl;
    //clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;
}

// Operator =
Geom24& Geom24::operator=(const Geom24& G)
{
    dim = G.get_dim();
    p = G.get_p();
    q = G.get_q();
    g2 = G.get_g2();
    nH = G.get_nH();
    nL = G.get_nL();
    nHL = G.get_nHL();
    dim_omega = G.get_dim_omega();

    // delete, reallocate and copy matrices
    delete [] mat;
    delete [] mom;
    delete [] omega;
    delete [] eps;
    mat = new arma::cx_mat [nHL];
    mom = new arma::cx_mat [nHL];
    eps = new int [nHL];
    omega = new arma::cx_mat [nHL];
    for(int i=0; i<nHL; i++)
    {
        mat[i] = G.get_mat(i);
        mom[i] = G.get_mat(i);
        eps[i] = G.get_eps(i);
        omega[i] = G.get_omega(i);
    }
    
    delete [] omega_table_4;
    int nHL4 = pow(nHL, 4);
    omega_table_4 = new cx_double [nHL4];
    for(int i=0; i<nHL4; i++)
        omega_table_4[i] = G.get_omega_table_4(i);
    
    //clog << "Geometry overwritten with the following parameters:" << endl;
    //clog << "(p, q, dim, g2) = (" << p << ", " << q << ", " << dim << ", " << g2 << ")" << endl;
    
    return *this;
}

// Destructor
Geom24::~Geom24()
{
    delete [] mat;
    delete [] mom;
    delete [] eps;
    delete [] omega;
    delete [] omega_table_4;
}


// Read parameters from istream
istream& Geom24::read_parameters(istream& in)
{
    if(in)
    {
        // read basic parameters
        in >> p >> q >> dim >> g2;
         
        // clear input stream state
        in.clear();
    }

    return in;
}


void Geom24::derived_parameters()
{
    int n = p+q;

    // create a type (p, q) clifford module
    Cliff C(p, q);
    vector<cx_mat> gamma = C.get_gamma();
    dim_omega = C.get_dim_gamma();
    
    vector<cx_mat> herm;
    vector<cx_mat> anti;

    for(int i=0; i<p; ++i)
        herm.push_back(gamma[i]);

    for(int i=0; i<q; ++i)
        anti.push_back(cx_double(0, 1)*gamma[p+i]);
    
    int  count = pow(2, n);
	// The outer for loop will run 2^n times (the number of all possible subsets).
	// Here variable i will act as a binary counter
	for (int i=0; i<count; i++)
	{
        vector<int> vec;
		// The inner for loop will run n times, As the maximum number of elements a set can have is n
		// This loop will generate a subset
		for (int j=0; j<n; j++)
		{
			// This if condition will check if jth bit in binary representation of i is set or not
			// if the value of (i & (1 << j)) is greater than 0, include arr[j] in the current subset
			// otherwise exclude arr[j]
			if ((i & (1 << j)) > 0)
                vec.push_back(j);
		}
        
        // Now calculate and push product if it has odd number of gammas
        int k = vec.size();
        if(k % 2 && k != 1)
        {
            vector<int>::const_iterator begin(vec.begin());
            vector<int>::const_iterator end(vec.end());
            cx_mat M = gamma.at(*begin);

            //cout << *begin << " ";
            for(vector<int>::const_iterator iter = vec.begin() + 1; iter != end; ++iter)
            {
                M *= gamma.at((*iter));
                //cout << *iter << " ";
            }

            if(M.is_hermitian())
            {
                herm.push_back(M);
                //cout << "herm" << endl;
            }
            else
            {
                anti.push_back(cx_double(0, 1)*M);
                //cout << "antiherm" << endl;
            }
        }
	}

    nH = herm.size();
    nL = anti.size();
    nHL = nH+nL;

    omega = new cx_mat [nHL];
    for(int i=0; i<nH; ++i)
        omega[i] = herm[i];
    for(int i=0; i<nL; ++i)
        omega[nH+i] = anti[i];

    init_omega_table_4();
}

ostream& operator<<(ostream& out, const Geom24& G)
{
    out << "Geometry (p, q, dim, g2) = (" << G.get_p() << ", " << G.get_q() << ", " << G.get_dim() << ", " << G.get_g2() << ") ";

    return out;
}

void Geom24::shuffle(gsl_rng* engine)
{
    for(int i=0; i<nHL; i++)
    {
        cx_mat temp(dim, dim);
        temp.imbue( [&engine](){ return cx_double(gsl_ran_gaussian(engine, 1.), gsl_ran_gaussian(engine, 1.)); } );

        mat[i] = (temp+temp.t())/(2*dim*dim);
    }
}

istream& Geom24::read_mat(istream& in)
{
    if(in)
    {
        // loop on matrices
        for(int i=0; i<nHL; ++i)
        {
            // loop on indices
            for(int j=0; j<dim; ++j)
            {
                for(int k=0; k<dim; ++k)
                {
                    double x, y;
                    in >> x >> y;
                    mat[i](j,k) = cx_double(x, y);
                }
            }
        }
         
        // clear input stream state
        in.clear();
    }

    return in;
}


void Geom24::reverse_mom()
{
    for(int i=0; i<nHL; ++i)
        mom[i] *= -1.;
}

void Geom24::init_omega_table_4()
{
    omega_table_4 = new cx_double [nHL*nHL*nHL*nHL];

    for(int i=0; i<nHL; ++i)
    {
        for(int j=0; j<nHL; ++j)
        {
            for(int k=0; k<nHL; ++k)
            {
                for(int l=0; l<nHL; ++l)
                    omega_table_4[l + nHL*(k + nHL*(j + nHL*i))] = trace(omega[i]*omega[j]*omega[k]*omega[l]);
            }
        }
    }
}

vector<int> base_conversion(int dec, const int& base, const int& max)
{
    vector<int> rem;

    while(dec)
    {
        rem.push_back(dec % base);
        dec /= base;
    }

    for(int i=rem.size(); i<max; i++)
        rem.push_back(0);

    reverse(rem.begin(), rem.end());

    return rem;
}


void Geom24::print_omega_table_4() const
{
    const int n = pow(nHL, 4);

    for(int i=0; i<n; ++i)
    {
        cx_double z = omega_table_4[i];
        if(z != cx_double(0., 0.))
        {
            int e = 1;
            vector<int> prod = base_conversion(i, nHL, 4);
            vector<int>::const_iterator end(prod.end());
            for(vector<int>::const_iterator iter = prod.begin(); iter != end; ++iter)
            {
                cout << (*iter) << " ";
                e *= eps[(*iter)];
            }
            cout << " " << omega_table_4[i] << e << endl;
        }
    }
}

void Geom24::adjust()
{
    // hermitianize
    for(int i=0; i<nHL; ++i)
        mat[i] = 0.5*(mat[i]+mat[i].t());

    // tracelessitize
    for(int i=nH; i<nHL; ++i)
    {
        double tr = trace(mat[i]).real()/dim;
        cx_mat idty(dim, dim, fill::eye);
        mat[i] = mat[i] - tr*idty;
    }
}
        





cx_mat Geom24::build_dirac() const
{
    // initialize dirac op to zero
    int dim_dirac = dim*dim*dim_omega;
    cx_mat dirac(dim_dirac, dim_dirac, fill::zeros);

    static cx_mat id(dim, dim, fill::eye);
    for(int i=0; i<nHL; ++i)
    {
        cx_mat bracket = kron(mat[i], id) + eps[i]*kron(id, mat[i].st());
        dirac += kron(omega[i], bracket);
    }

    return dirac;
}

double Geom24::calculate_S_from_dirac() const
{
    cx_mat dirac = build_dirac();
    cx_mat dirac2 = dirac*dirac;
    double trdirac2 = trace(dirac2).real();
    double trdirac4 = trace(dirac2*dirac2).real();
    return g2*trdirac2 + trdirac4;
}

double Geom24::dirac2() const
{
    double res = 0.;
    for(int i=0; i<nHL; ++i)
    {
        double tr2 = trace(mat[i]*mat[i]).real();
        double tr1 = trace(mat[i]).real();

        res += (dim*tr2 + eps[i]*tr1*tr1);
    }

    return 2.*dim_omega*res;
}




double Geom24::calculate_S() const
{
    return g2*dirac2() + dirac4();
    //return dirac2();
}


double Geom24::compute_A4(const int& i1, const int& i2, const int& i3, const int& i4) const
{
    // epsilon factor
    int e = eps[i1]*eps[i2]*eps[i3]*eps[i4];

    // if e=-1, then [1+*e] becomes 2i*imag
    // and the clifford part is guaranteed to
    // be pure imaginary
    if(e<0)
    {
        // clifford product
        double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].imag(); 

        if(cliff != 0.)
        {
            // base matrix products
            cx_mat M1M2 = mat[i1]*mat[i2];
            cx_mat M1M3 = mat[i1]*mat[i3];
            cx_mat M1M4 = mat[i1]*mat[i4];
            cx_mat M2M3 = mat[i2]*mat[i3];
            cx_mat M2M4 = mat[i2]*mat[i4];
            cx_mat M3M4 = mat[i3]*mat[i4];

            // traces
            double tr1234 = trace(M1M2*M3M4).imag();
            double tr234 = trace(M2M3*mat[i4]).imag();
            double tr134 = trace(M1M3*mat[i4]).imag();
            double tr124 = trace(M1M2*mat[i4]).imag();
            double tr123 = trace(M1M2*mat[i3]).imag();
            double tr1 = trace(mat[i1]).real();
            double tr2 = trace(mat[i2]).real();
            double tr3 = trace(mat[i3]).real();
            double tr4 = trace(mat[i4]).real();
            
            // compute sum
            double res = dim*tr1234;
            res += eps[i1]*tr1*tr234;
            res += eps[i2]*tr2*tr134;
            res += eps[i3]*tr3*tr124;
            res += eps[i4]*tr4*tr123;

            return -2*cliff*res;
            // NOTE: this minus here comes from the i in cliff
            // and the i coming from 2i*imag
        }
        else
            return 0.;
    }
    else
    {
        // clifford product
        double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].real(); 

        if(cliff != 0.)
        {
            // base matrix products
            cx_mat M1M2 = mat[i1]*mat[i2];
            cx_mat M1M3 = mat[i1]*mat[i3];
            cx_mat M1M4 = mat[i1]*mat[i4];
            cx_mat M2M3 = mat[i2]*mat[i3];
            cx_mat M2M4 = mat[i2]*mat[i4];
            cx_mat M3M4 = mat[i3]*mat[i4];

            // traces
            double tr1234 = trace(M1M2*M3M4).real();
            double tr234 = trace(M2M3*mat[i4]).real();
            double tr134 = trace(M1M3*mat[i4]).real();
            double tr124 = trace(M1M2*mat[i4]).real();
            double tr123 = trace(M1M2*mat[i3]).real();
            double tr12 = trace(M1M2).real();
            double tr34 = trace(M3M4).real();
            double tr13 = trace(M1M3).real();
            double tr24 = trace(M2M4).real();
            double tr14 = trace(M1M4).real();
            double tr23 = trace(M2M3).real();
            double tr1 = trace(mat[i1]).real();
            double tr2 = trace(mat[i2]).real();
            double tr3 = trace(mat[i3]).real();
            double tr4 = trace(mat[i4]).real();


            double res = dim*tr1234;
            res += eps[i1]*tr1*tr234;
            res += eps[i2]*tr2*tr134;
            res += eps[i3]*tr3*tr124;
            res += eps[i4]*tr4*tr123;
            res += eps[i1]*eps[i2]*tr12*tr34;
            res += eps[i1]*eps[i3]*tr13*tr24;
            res += eps[i1]*eps[i4]*tr14*tr23;

            double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))].real(); 

            return 2*cliff*res;
        }
        else
            return 0.;
    }
}

double Geom24::compute_A2(const int& i1, const int& i2) const
{
    // clifford product
    double cliff = omega_table_4[i2 + nHL*(i1 + nHL*(i2 + nHL*i1))].real();

    // base matrix products
    cx_mat M1M1 = mat[i1]*mat[i1];
    cx_mat M2M2 = mat[i2]*mat[i2];
    cx_mat M1M2 = mat[i1]*mat[i2];

    // traces
    double tr1122 = trace(M1M1*M2M2).real();
    double tr1212 = trace(M1M2*M1M2).real();
    double tr122 = trace(M1M2*mat[i2]).real();
    double tr112 = trace(M1M1*mat[i2]).real();
    double tr11 = trace(M1M1).real();
    double tr22 = trace(M2M2).real();
    double tr12 = trace(M1M2).real();
    double tr1 = trace(mat[i1]).real();
    double tr2 = trace(mat[i2]).real();
    
    
    if(cliff < 0)
    {
        // compute sum
        double res = dim*(2*tr1122 - tr1212);
        res += 2*eps[i1]*tr1*tr122;
        res += 2*eps[i2]*tr2*tr112;
        res += tr11*tr22;
        res += 2*eps[i1]*eps[i2]*tr12*tr12;

        return 2*dim_omega*res;
    }
    else
    {
        // compute sum
        double res = dim*(2*tr1122 + tr1212);
        res += 6*eps[i1]*tr1*tr122;
        res += 6*eps[i2]*tr2*tr112;
        res += 3*tr11*tr22;
        res += 6*eps[i1]*eps[i2]*tr12*tr12;

        return 2*dim_omega*res;
    }
}

double Geom24::compute_A(const int& i) const
{
    // base matrix products
    cx_mat M2 = mat[i]*mat[i];
    cx_mat M3 = M2*mat[i];

    // traces
    double tr1 = trace(mat[i]).real();
    double tr2 = trace(M2).real();
    double tr3 = trace(M3).real();
    double tr4 = trace(M3*mat[i]).real();

    double res = dim*tr4;
    res += 4*eps[i]*tr1*tr3;
    res += 3*tr2*tr2;

    return 2*dim_omega*res;
}

double Geom24::dirac4() const
{
    double res = 0.;

    // four distinct indices
    for(int i=0; i<nHL; ++i)
    {
        for(int j=i+1; j<nHL; ++j)
        {
            for(int k=j+1; k<nHL; ++k)
            {
                for(int l=k+1; l<nHL; ++l)
                    res += 8*(compute_A4(i,j,k,l)+compute_A4(i,j,l,k)+compute_A4(i,k,j,l));
            }
        }
    }

    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        for(int j=i+1; j<nHL; ++j)
            res += 2*compute_A2(i,j);
    }

    // all indices equal
    for(int i=0; i<nHL; ++i)
        res += compute_A(i);

    return res;
}


void Geom24::sample_mom(gsl_rng* engine)
{
    for(int i=0; i<nHL; ++i)
    {
        cx_mat temp(dim, dim);
        temp.imbue( [&engine](){ return cx_double(gsl_ran_gaussian(engine, 1.), gsl_ran_gaussian(engine, 1.)); } );

        mom[i] = (temp+temp.t())/2.;
    }
}

double Geom24::calculate_K() const
{
    double res = 0;

    for(int i=0; i<nHL; ++i)
        res += trace(mom[i]*mom[i]).real();

    return res/2;
}

double Geom24::calculate_H() const
{
    return calculate_S() + calculate_K();
}

cx_mat Geom24::compute_B4(const int& k, const int& i2, const int& i3, const int& i4, const double& cliff, const bool& neg) const
{
    if(neg)
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        double tr234 = trace(M2M3M4).imag();
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cx_double iu(0,1);
        cx_mat res(dim ,dim, fill::eye);
        res *= -2*eps[k]*tr234;
        res += double(dim)*iu*(M2M3M4 - M2M3M4.t());
        res += eps[i2]*tr2*iu*(M3M4 - M3M4.t());
        res += eps[i3]*tr3*iu*(M2M4 - M2M4.t());
        res += eps[i4]*tr4*iu*(M2M3 - M2M3.t());

        return cliff*res.st();
    }
    else
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        double tr234 = trace(M2M3M4).real();
        double tr23 = trace(M2M3).real();
        double tr24 = trace(M2M4).real();
        double tr34 = trace(M3M4).real();
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cx_mat res(dim ,dim, fill::eye);
        res *= 2*eps[k]*tr234;
        res += dim*(M2M3M4 + M2M3M4.t());
        res += eps[i2]*tr2*(M3M4 + M3M4.t());
        res += eps[i3]*tr3*(M2M4 + M2M4.t());
        res += eps[i4]*tr4*(M2M3 + M2M3.t());
        res += 2*eps[k]*eps[i2]*tr34*mat[i2];
        res += 2*eps[k]*eps[i3]*tr24*mat[i3];
        res += 2*eps[k]*eps[i4]*tr23*mat[i4];

        return cliff*res.st();
    }
}

void Geom24::compute_B4_cout(const int& k, const int& i2, const int& i3, const int& i4, const double& cliff, const bool& neg) const
{
    if(neg)
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        double tr234 = trace(M2M3M4).imag();
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cout << -2*cliff*eps[k]*tr234 << endl;
        cout << -double(dim)*cliff*2*M2M3M4(0,0).imag() << endl;
        cout << -eps[i2]*tr2*cliff*2*M3M4(0,0).imag() << endl;
        cout << -eps[i3]*tr3*cliff*2*M2M4(0,0).imag() << endl;
        cout << -eps[i4]*tr4*cliff*2*M2M3(0,0).imag() << endl;
    }
    else
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        double tr234 = trace(M2M3M4).real();
        double tr23 = trace(M2M3).real();
        double tr24 = trace(M2M4).real();
        double tr34 = trace(M3M4).real();
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cout << 2*cliff*eps[k]*tr234 << endl;
        cout << dim*2*cliff*M2M3M4(0,0).real() << endl;
        cout << eps[i2]*cliff*tr2*2*M3M4(0,0).real() << endl;
        cout << eps[i3]*cliff*tr3*2*M2M4(0,0).real() << endl;
        cout << eps[i4]*cliff*tr4*2*M2M3(0,0).real() << endl;
        cout << 2*eps[k]*eps[i2]*cliff*tr34*mat[i2](0,0).real() << endl;
        cout << 2*eps[k]*eps[i3]*cliff*tr24*mat[i3](0,0).real() << endl;
        cout << 2*eps[k]*eps[i4]*cliff*tr23*mat[i4](0,0).real() << endl;
    }
}

cx_mat Geom24::compute_B4_bruteforce(const int& k, const int& i2, const int& i3, const int& i4, const cx_double& cliff, const int& e) const
{
    // base matrix products
    cx_mat M2M3 = mat[i2]*mat[i3];
    cx_mat M2M4 = mat[i2]*mat[i4];
    cx_mat M3M4 = mat[i3]*mat[i4];
    cx_mat M2M3M4 = M2M3*mat[i4];

    // traces
    cx_double tr234 = trace(M2M3M4);
    double tr23 = trace(M2M3).real();
    double tr24 = trace(M2M4).real();
    double tr34 = trace(M3M4).real();
    double tr2 = trace(mat[i2]).real();
    double tr3 = trace(mat[i3]).real();
    double tr4 = trace(mat[i4]).real();

    // compute sum
    cx_mat res(dim ,dim, fill::eye);
    res *= eps[k];
    res *= tr234 + double(e)*conj(tr234);
    res += dim*(M2M3M4 + e*M2M3M4.t());
    res += eps[i2]*tr2*(M3M4 + e*M3M4.t());
    res += eps[i3]*tr3*(M2M4 + e*M2M4.t());
    res += eps[i4]*tr4*(M2M3 + e*M2M3.t());
    res += (1+e)*eps[k]*eps[i2]*tr34*mat[i2];
    res += (1+e)*eps[k]*eps[i3]*tr24*mat[i3];
    res += (1+e)*eps[k]*eps[i4]*tr23*mat[i4];

    return cliff*res.st();
}

cx_mat Geom24::compute_B4_explicit(const int& k, const int& i2, const int& i3, const int& i4, const bool& neg) const
{
    if(neg)
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        cx_double tr234 = trace(M2M3M4);
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cx_mat idty(dim ,dim, fill::eye);
        cx_mat res = double(eps[k])*tr234*idty;
        res += -double(eps[k])*conj(tr234)*idty;
        res += dim*(M2M3M4 - M2M3M4.t());
        res += eps[i2]*tr2*(M3M4 - M3M4.t());
        res += eps[i3]*tr3*(M2M4 - M2M4.t());
        res += eps[i4]*tr4*(M2M3 - M2M3.t());

        return res.st();
    }
    else
    {
        // base matrix products
        cx_mat M2M3 = mat[i2]*mat[i3];
        cx_mat M2M4 = mat[i2]*mat[i4];
        cx_mat M3M4 = mat[i3]*mat[i4];
        cx_mat M2M3M4 = M2M3*mat[i4];

        // traces
        cx_double tr234 = trace(M2M3M4);
        double tr23 = trace(M2M3).real();
        double tr24 = trace(M2M4).real();
        double tr34 = trace(M3M4).real();
        double tr2 = trace(mat[i2]).real();
        double tr3 = trace(mat[i3]).real();
        double tr4 = trace(mat[i4]).real();

        // compute sum
        cx_mat idty(dim ,dim, fill::eye);
        cx_mat res = double(eps[k])*tr234*idty;
        res += double(eps[k])*conj(tr234)*idty;
        res += dim*(M2M3M4 + M2M3M4.t());
        res += eps[i2]*tr2*(M3M4 + M3M4.t());
        res += eps[i3]*tr3*(M2M4 + M2M4.t());
        res += eps[i4]*tr4*(M2M3 + M2M3.t());
        res += 2*eps[k]*eps[i2]*tr34*mat[i2];
        res += 2*eps[k]*eps[i3]*tr24*mat[i3];
        res += 2*eps[k]*eps[i4]*tr23*mat[i4];

        return res.st();
    }
}


cx_mat Geom24::compute_B2(const int& k, const int& i) const
{
    // clifford product
    double cliff = omega_table_4[i + nHL*(k + nHL*(i + nHL*k))].real();

    // base matrix products
    cx_mat MiMk = mat[i]*mat[k];
    cx_mat MiMi = mat[i]*mat[i];
    cx_mat MiMiMk = mat[i]*MiMk;
    cx_mat MiMkMi = MiMk*mat[i];

    // traces
    double triki = trace(MiMkMi).real();
    double trik = trace(MiMk).real();
    double trii = trace(MiMi).real();
    double tri = trace(mat[i]).real();
    double trk = trace(mat[k]).real();
    
    
    if(cliff < 0)
    {
        // compute sum
        cx_mat res(dim, dim, fill::eye);
        res *= eps[k]*triki;
        res += dim*(MiMiMk + MiMiMk.t() - MiMkMi);
        res += eps[i]*tri*(MiMk + MiMk.t());
        res += 2*eps[k]*eps[i]*trik*mat[i];
        res += eps[k]*trk*MiMi;
        res += trii*mat[k];

        return 2*dim_omega*res.st();
    }
    else
    {
        // compute sum
        cx_mat res(dim, dim, fill::eye);
        res *= 3*eps[k]*triki;
        res += dim*(MiMiMk + MiMiMk.t() + MiMkMi);
        res += 3*eps[i]*tri*(MiMk + MiMk.t());
        res += 6*eps[k]*eps[i]*trik*mat[i];
        res += 3*eps[k]*trk*MiMi;
        res += 3*trii*mat[k];

        return 2*dim_omega*res.st();
    }
}

void Geom24::compute_B2_cout(const int& k, const int& i) const
{
    // clifford product
    double cliff = omega_table_4[i + nHL*(k + nHL*(i + nHL*k))].real();

    // base matrix products
    cx_mat MiMk = mat[i]*mat[k];
    cx_mat MiMi = mat[i]*mat[i];
    cx_mat MiMiMk = mat[i]*MiMk;
    cx_mat MiMkMi = MiMk*mat[i];

    // traces
    double triki = trace(MiMkMi).real();
    double trik = trace(MiMk).real();
    double trii = trace(MiMi).real();
    double tri = trace(mat[i]).real();
    double trk = trace(mat[k]).real();
    
    
    if(cliff < 0)
    {
        // compute sum
        cout << eps[k]*2*dim_omega*triki << endl;
        cout << dim*2*dim_omega*(2*MiMiMk(0,0).real() - MiMkMi(0,0).real()) << endl;
        cout << eps[i]*tri*2*dim_omega*2*MiMk(0,0).real() << endl;
        cout << 2*eps[k]*eps[i]*trik*2*dim_omega*mat[i](0,0).real() << endl;
        cout << eps[k]*trk*2*dim_omega*MiMi(0,0).real() << endl;
        cout << trii*2*dim_omega*mat[k](0,0).real() << endl;
    }
    else
    {
        // compute sum
        cout << 3*eps[k]*2*dim_omega*triki << endl;
        cout << dim*2*dim_omega*(2*MiMiMk(0,0).real() + MiMkMi(0,0).real()) << endl;
        cout << 3*eps[i]*2*dim_omega*tri*2*MiMk(0,0).real() << endl;
        cout << 6*eps[k]*eps[i]*2*dim_omega*trik*mat[i](0,0).real() << endl;
        cout << 3*eps[k]*2*dim_omega*trk*MiMi(0,0).real() << endl;
        cout << 3*2*dim_omega*trii*mat[k](0,0).real() << endl;
    }
}

cx_mat Geom24::compute_B2_iki_explicit(const int& k, const int& i) const
{
    // base matrix products
    cx_mat MiMk = mat[i]*mat[k];
    cx_mat MiMi = mat[i]*mat[i];
    cx_mat MiMkMi = MiMk*mat[i];

    // traces
    double triki = trace(MiMkMi).real();
    double trik = trace(MiMk).real();
    double trii = trace(MiMi).real();
    double tri = trace(mat[i]).real();
    double trk = trace(mat[k]).real();
    
    
    // compute sum
    cx_mat res(dim, dim, fill::eye);
    res *= eps[k]*triki;
    res += dim*MiMkMi;
    res += eps[i]*tri*(MiMk + MiMk.t());
    res += 2*eps[k]*eps[i]*trik*mat[i];
    res += eps[k]*trk*MiMi;
    res += trii*mat[k];

    return 2*res.st();
}

cx_mat Geom24::compute_B2_iik_explicit(const int& k, const int& i) const
{
    // base matrix products
    cx_mat MiMk = mat[i]*mat[k];
    cx_mat MiMi = mat[i]*mat[i];
    cx_mat MiMiMk = mat[i]*MiMk;
    cx_mat MiMkMi = MiMk*mat[i];

    // traces
    double triik = trace(MiMiMk).real();
    double trik = trace(MiMk).real();
    double trii = trace(MiMi).real();
    double tri = trace(mat[i]).real();
    double trk = trace(mat[k]).real();
    
    
    // compute sum
    cx_mat res(dim, dim, fill::eye);
    res *= eps[k]*triik;
    res += 0.5*dim*(MiMiMk + MiMiMk.t());
    res += eps[i]*tri*(MiMk + MiMk.t());
    res += 2*eps[k]*eps[i]*trik*mat[i];
    res += eps[k]*trk*MiMi;
    res += trii*mat[k];

    return 2*res.st();
}


cx_mat Geom24::compute_B(const int& k) const
{
    // base matrix products
    cx_mat M2 = mat[k]*mat[k];
    cx_mat M3 = mat[k]*M2;

    // traces
    double tr3 = trace(M3).real();
    double tr2 = trace(M2).real();
    double tr1 = trace(mat[k]).real();

    cx_mat res(dim, dim, fill::eye);
    res *= eps[k]*tr3;
    res += dim*M3;
    res += 3*tr2*mat[k];
    res += 3*eps[k]*tr1*M2;

    return 2*dim_omega*res.st();
}

void Geom24::compute_B_cout(const int& k) const
{
    // base matrix products
    cx_mat M2 = mat[k]*mat[k];
    cx_mat M3 = mat[k]*M2;

    // traces
    double tr3 = trace(M3).real();
    double tr2 = trace(M2).real();
    double tr1 = trace(mat[k]).real();

    cout << eps[k]*2*dim_omega*tr3 << endl;
    cout << dim*2*dim_omega*M3(0,0).real() << endl;
    cout << 3*2*dim_omega*tr2*mat[k](0,0).real() << endl;
    cout << 3*eps[k]*2*dim_omega*tr1*M2(0,0).real() << endl;
}

cx_mat Geom24::der_dirac4(const int& k, const bool& herm) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        if(i1 != k)
        {
            for(int i2=i1+1; i2<nHL; ++i2)
            {
                if(i2 != k)
                {
                    for(int i3=i2+1; i3<nHL; ++i3)
                    {
                        if(i3 != k)
                        {
                            // epsilon factor
                            int e = eps[k]*eps[i1]*eps[i2]*eps[i3];

                            if(e<0)
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].imag(); 

                                if(cliff1 != 0.)
                                {
                                    cx_mat temp = compute_B4(k,i1,i2,i3, cliff1, true) + compute_B4(k,i1,i3,i2, cliff2, true) + compute_B4(k,i2,i1,i3, cliff3, true);
                                    res += temp + temp.t();
                                }
                            }
                            else
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].real(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].real(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].real(); 

                                if(cliff1 != 0.)
                                {
                                    cx_mat temp = compute_B4(k,i1,i2,i3, cliff1, false) + compute_B4(k,i1,i3,i2, cliff2, false) + compute_B4(k,i2,i1,i3, cliff3, false);
                                    res += temp + temp.t();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        if(i != k)
            res += compute_B2(k,i);
    }

    // all indices equal
    res += compute_B(k);


    if(herm)
        return 2*(res+res.t());
    else
        return 4*res;

}

void Geom24::der_dirac4_cout(const int& k) const
{
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        if(i1 != k)
        {
            for(int i2=i1+1; i2<nHL; ++i2)
            {
                if(i2 != k)
                {
                    for(int i3=i2+1; i3<nHL; ++i3)
                    {
                        if(i3 != k)
                        {
                            // epsilon factor
                            int e = eps[k]*eps[i1]*eps[i2]*eps[i3];

                            if(e<0)
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].imag(); 

                                if(cliff1 != 0.)
                                {
                                    compute_B4_cout(k,i1,i2,i3, cliff1, true);
                                    compute_B4_cout(k,i1,i3,i2, cliff2, true);
                                    compute_B4_cout(k,i2,i1,i3, cliff3, true);
                                }
                            }
                            else
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].real(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].real(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].real(); 

                                if(cliff1 != 0.)
                                {
                                    compute_B4_cout(k,i1,i2,i3, cliff1, false);
                                    compute_B4_cout(k,i1,i3,i2, cliff2, false);
                                    compute_B4_cout(k,i2,i1,i3, cliff3, false);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        if(i != k)
            compute_B2_cout(k,i);
    }

    // all indices equal
    compute_B_cout(k);
}

cx_mat Geom24::der_dirac4_bruteforce(const int& k, const bool& herm) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        for(int i2=0; i2<nHL; ++i2)
        {
            for(int i3=0; i3<nHL; ++i3)
            {
                for(int i4=0; i4<nHL; ++i4)
                {
                    // epsilon factor
                    int e = eps[i1]*eps[i2]*eps[i3]*eps[i4];

                    // clifford product
                    cx_double cliff = omega_table_4[i4 + nHL*(i3 + nHL*(i2 + nHL*i1))]; 

                    // compute derivative
                    if(k==i1)
                        res += compute_B4_bruteforce(i1,i2,i3,i4, cliff, e);
                    if(k==i2)
                        res += compute_B4_bruteforce(i2,i3,i4,i1, cliff, e);
                    if(k==i3)
                        res += compute_B4_bruteforce(i3,i4,i1,i2, cliff, e);
                    if(k==i4)
                        res += compute_B4_bruteforce(i4,i1,i2,i3, cliff, e);
                }
            }
        }
    }

    if(herm)
        return 0.5*(res+res.t());
    else
        return res;
}


cx_mat Geom24::der_dirac4_explicit(const int& k, const bool& herm) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        if(i1 != k)
        {
            for(int i2=i1+1; i2<nHL; ++i2)
            {
                if(i2 != k)
                {
                    for(int i3=i2+1; i3<nHL; ++i3)
                    {
                        if(i3 != k)
                        {
                            // epsilon factor
                            int e = eps[k]*eps[i1]*eps[i2]*eps[i3];
                            
                            // clifford product
                            cx_double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))]; 
                            cx_double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))]; 
                            cx_double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))]; 

                            if(e<0)
                            {
                                cx_mat temp = cliff1*compute_B4_explicit(k,i1,i2,i3, true) + cliff2*compute_B4_explicit(k,i1,i3,i2, true) + cliff3*compute_B4_explicit(k,i2,i1,i3, true);
                                res += temp + temp.t();
                            }
                            else if(e>0)
                            {
                                cx_mat temp = cliff1*compute_B4_explicit(k,i1,i2,i3, false) + cliff2*compute_B4_explicit(k,i1,i3,i2, false) + cliff3*compute_B4_explicit(k,i2,i1,i3, false);
                                res += temp + temp.t();
                            }
                        }
                    }
                }
            }
        }
    }

    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        if(i != k)
        {
            cx_double cliff = omega_table_4[i + nHL*(k + nHL*(i + nHL*k))];
            res += 2*dim_omega*compute_B2_iik_explicit(k,i);
            res += cliff*compute_B2_iki_explicit(k,i);
        }
    }

    // all indices equal
    res += compute_B(k);


    if(herm)
        return 2*(res+res.t());
    else
        return 4*res;
}

cx_mat Geom24::debug_4different(const int& k) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        if(i1 != k)
        {
            for(int i2=i1+1; i2<nHL; ++i2)
            {
                if(i2 != k)
                {
                    for(int i3=i2+1; i3<nHL; ++i3)
                    {
                        if(i3 != k)
                        {
                            // epsilon factor
                            int e = eps[k]*eps[i1]*eps[i2]*eps[i3];

                            if(e<0)
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].imag(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].imag(); 

                                if(cliff1 != 0.)
                                {
                                    cx_mat temp = compute_B4(k,i1,i2,i3, cliff1, true) + compute_B4(k,i1,i3,i2, cliff2, true) + compute_B4(k,i2,i1,i3, cliff3, true);
                                    res += temp + temp.t();
                                }
                            }
                            else
                            {
                                // clifford product
                                double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))].real(); 
                                double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))].real(); 
                                double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))].real(); 

                                if(cliff1 != 0.)
                                {
                                    cx_mat temp = compute_B4(k,i1,i2,i3, cliff1, false) + compute_B4(k,i1,i3,i2, cliff2, false) + compute_B4(k,i2,i1,i3, cliff3, false);
                                    res += temp + temp.t();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 4*res;

}

cx_mat Geom24::debug_4different_explicit(const int& k) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // four distinct indices
    for(int i1=0; i1<nHL; ++i1)
    {
        if(i1 != k)
        {
            for(int i2=i1+1; i2<nHL; ++i2)
            {
                if(i2 != k)
                {
                    for(int i3=i2+1; i3<nHL; ++i3)
                    {
                        if(i3 != k)
                        {
                            // epsilon factor
                            int e = eps[k]*eps[i1]*eps[i2]*eps[i3];
                            
                            // clifford product
                            cx_double cliff1 = omega_table_4[i3 + nHL*(i2 + nHL*(i1 + nHL*k))]; 
                            cx_double cliff2 = omega_table_4[i2 + nHL*(i3 + nHL*(i1 + nHL*k))]; 
                            cx_double cliff3 = omega_table_4[i3 + nHL*(i1 + nHL*(i2 + nHL*k))]; 

                            if(e<0)
                            {
                                cx_mat temp = cliff1*compute_B4_explicit(k,i1,i2,i3, true) + cliff2*compute_B4_explicit(k,i1,i3,i2, true) + cliff3*compute_B4_explicit(k,i2,i1,i3, true);
                                res += temp + temp.t();
                            }
                            else if(e>0)
                            {
                                cx_mat temp = cliff1*compute_B4_explicit(k,i1,i2,i3, false) + cliff2*compute_B4_explicit(k,i1,i3,i2, false) + cliff3*compute_B4_explicit(k,i2,i1,i3, false);
                                res += temp + temp.t();
                            }
                        }
                    }
                }
            }
        }
    }

    return 4*res;
}

cx_mat Geom24::debug_2equal_explicit(const int& k) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        if(i != k)
        {
            cx_double cliff = omega_table_4[i + nHL*(k + nHL*(i + nHL*k))];
            res += 2*dim_omega*compute_B2_iik_explicit(k,i);
            res += cliff*compute_B2_iki_explicit(k,i);
        }
    }

    return res;
}

cx_mat Geom24::debug_2equal(const int& k) const
{
    cx_mat res(dim, dim, fill::zeros);
    
    // two distinct pairs of equal indices
    for(int i=0; i<nHL; ++i)
    {
        if(i != k)
            res += compute_B2(k,i);
    }

    return res;
}

cx_mat Geom24::der_dirac2(const int& k) const
{
    cx_mat res(dim, dim, fill::eye);

    res *= eps[k]*trace(mat[k]).real();
    res += dim*mat[k].st();

    return 4*dim_omega*res;
}

void Geom24::der_dirac2_cout(const int& k) const
{
    cout << g2*4*dim_omega*eps[k]*trace(mat[k]).real() << endl;
    cout << g2*4*dim_omega*dim*mat[k](0,0).real() << endl;
}

cx_mat Geom24::der_dirac24(const int& k, const bool& herm) const
{
    return g2*der_dirac2(k) + der_dirac4(k, herm);
}


// old function for trD4

double Geom24::dirac4_old() const
{
    double res = 0.;
    int* i = new int [4];


    for(i[3]=0; i[3]<nHL; i[3]++)
    {
        for(i[2]=0; i[2]<nHL; i[2]++)
        {
            for(i[1]=0; i[1]<=i[2]; i[1]++)
            {
                for(i[0]=0; i[0]<=i[3]; i[0]++)
                {
                    if(i[0]==i[1] && i[1]==i[2] && i[2]==i[3])
                    {
                        // compute matrix products
                        cx_mat M0M0 = mat[i[0]]*mat[i[0]];
                        cx_mat M0M0M0 = M0M0*mat[i[0]];
                        cx_mat M0M0M0M0 = M0M0M0*mat[i[0]];

                        // compute traces
                        double trM0 = trace(mat[i[0]]).real();
                        double trM0M0 = trace(M0M0).real();
                        double trM0M0M0 = trace(M0M0M0).real();
                        double trM0M0M0M0 = trace(M0M0M0M0).real();
                        
                        // add to total
                        double temp = 0;

                        // tr4 term
                        temp += dim*2.*trM0M0M0M0;

                        // tr3tr1 term
                        temp += 8.*eps[i[0]]*trM0M0M0*trM0;

                        // tr2tr2 term
                        temp += 6.*trM0M0*trM0M0;
                        
                        res += dim_omega*temp;
                    }


                    
                    else if(i[0] == i[3] && i[1] == i[2])
                    {
                        // alloc matrix products (trace needed)
                        cx_mat M0M0 = mat[i[0]]*mat[i[0]];
                        cx_mat M1M1 = mat[i[1]]*mat[i[1]];
                        cx_mat M0M1 = mat[i[0]]*mat[i[1]];
                        cx_mat M1M0M0M1 = mat[i[1]]*M0M0*mat[i[1]];
                        cx_mat M0M1M0 = M0M1*mat[i[0]];
                        cx_mat M1M0M1 = mat[i[1]]*M0M1;
                            
                        // compute traces
                        double trM0 = trace(mat[i[0]]).real();
                        double trM1 = trace(mat[i[1]]).real();
                        double trM0M0 = trace(M0M0).real();
                        double trM1M1 = trace(M1M1).real();
                        double trM0M1 = trace(M0M1).real();
                        double trM0M1M0 = trace(M0M1M0).real();
                        double trM1M0M1 = trace(M1M0M1).real();
                        double trM1M0M0M1 = trace(M1M0M0M1).real();
                        
                        // add to total
                        double temp = 0;

                        // tr4 term
                        temp += dim*2.*trM1M0M0M1;

                        // tr3tr1 term
                        temp += 4.*eps[i[0]]*trM1M0M1*trM0;
                        temp += 4.*eps[i[1]]*trM0M1M0*trM1;

                        // tr2tr2 term
                        temp += 4.*eps[i[0]]*eps[i[1]]*trM0M1*trM0M1;
                        temp += 2.*trM0M0*trM1M1;
                        
                        res += dim_omega*temp;
                    }
                    

                    else
                    {
                        cx_double cliff = omega_table_4[i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                        
                        if(cliff.real() != 0. || cliff.imag() != 0.)
                        {
                            // compute matrix products
                            cx_mat M0M1 = mat[i[0]]*mat[i[1]];
                            cx_mat M0M2 = mat[i[0]]*mat[i[2]];
                            cx_mat M0M3 = mat[i[0]]*mat[i[3]];
                            cx_mat M1M2 = mat[i[1]]*mat[i[2]];
                            cx_mat M1M3 = mat[i[1]]*mat[i[3]];
                            cx_mat M2M3 = mat[i[2]]*mat[i[3]];
                            cx_mat M0M1M2 = M0M1*mat[i[2]];
                            cx_mat M0M1M3 = M0M1*mat[i[3]];
                            cx_mat M0M2M3 = M0M2*mat[i[3]];
                            cx_mat M1M2M3 = mat[i[1]]*M2M3;
                            cx_mat M0M1M2M3 = mat[i[0]]*M1M2M3;

                            // compute traces
                            cx_double trM0M1M2M3 = trace(M0M1M2M3);
                            cx_double trM0M1M2 = trace(M0M1M2);
                            cx_double trM0M1M3 = trace(M0M1M3);
                            cx_double trM0M2M3 = trace(M0M2M3);
                            cx_double trM1M2M3 = trace(M1M2M3);
                            double trM0M1 = trace(M0M1).real();
                            double trM0M2 = trace(M0M2).real();
                            double trM0M3 = trace(M0M3).real();
                            double trM1M2 = trace(M1M2).real();
                            double trM1M3 = trace(M1M3).real();
                            double trM2M3 = trace(M2M3).real();
                            double trM0 = trace(mat[i[0]]).real();
                            double trM1 = trace(mat[i[1]]).real();
                            double trM2 = trace(mat[i[2]]).real();
                            double trM3 = trace(mat[i[3]]).real();
                            

                            // add to total

                            // tr4 terms
                            cx_double T1 = trM0M1M2M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM0M1M2M3);
                            res += 2.*dim*(cliff*T1).real();

                            // tr3tr1 terms

                            cx_double T3 = (double)(eps[i[3]])*trM0M1M2 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]])*conj(trM0M1M2);
                            T3 = trM3*T3;

                            cx_double T4 = (double)(eps[i[2]])*trM0M1M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[3]])*conj(trM0M1M3);
                            T3 += trM2*T4;

                            cx_double T5 = (double)(eps[i[1]])*trM0M2M3 + (double)(eps[i[0]]*eps[i[2]]*eps[i[3]])*conj(trM0M2M3);
                            T3 += trM1*T5;

                            cx_double T6 = (double)(eps[i[0]])*trM1M2M3 + (double)(eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM1M2M3);
                            T3 += trM0*T6;

                            res += 2.*(cliff*T3).real();

                            // tr2tr2 terms
                            double T7 = trM0M1*trM2M3*(eps[i[0]]*eps[i[1]] + eps[i[2]]*eps[i[3]]);
                            double T8 = trM0M2*trM1M3*(eps[i[0]]*eps[i[2]] + eps[i[1]]*eps[i[3]]);
                            double T9 = trM0M3*trM1M2*(eps[i[0]]*eps[i[3]] + eps[i[1]]*eps[i[2]]);

                            res += 2.*cliff.real()*(T7+T8+T9);
                        }
                    }
                }
            }
        }
    }
    
    for(i[3]=0; i[3]<nHL; i[3]++)
    {
        for(i[1]=0; i[1]<nHL; i[1]++)
        {
            for(i[2]=0; i[2]<i[1]; i[2]++)
            {
                for(i[0]=0; i[0]<i[3]; i[0]++)
                {
                    cx_double cliff = omega_table_4[i[3] + nHL*(i[2] + nHL*(i[1] + nHL*i[0]))];
                    
                    if(cliff.real() != 0. || cliff.imag() != 0.)
                    {
                        // compute matrix products
                        cx_mat M0M1 = mat[i[0]]*mat[i[1]];
                        cx_mat M0M2 = mat[i[0]]*mat[i[2]];
                        cx_mat M0M3 = mat[i[0]]*mat[i[3]];
                        cx_mat M1M2 = mat[i[1]]*mat[i[2]];
                        cx_mat M1M3 = mat[i[1]]*mat[i[3]];
                        cx_mat M2M3 = mat[i[2]]*mat[i[3]];
                        cx_mat M0M1M2 = M0M1*mat[i[2]];
                        cx_mat M0M1M3 = M0M1*mat[i[3]];
                        cx_mat M0M2M3 = M0M2*mat[i[3]];
                        cx_mat M1M2M3 = M1M2*mat[i[3]];
                        cx_mat M0M1M2M3 = mat[i[0]]*M1M2M3;

                        // compute traces
                        cx_double trM0M1M2M3 = trace(M0M1M2M3);
                        cx_double trM0M1M2 = trace(M0M1M2);
                        cx_double trM0M1M3 = trace(M0M1M3);
                        cx_double trM0M2M3 = trace(M0M2M3);
                        cx_double trM1M2M3 = trace(M1M2M3);
                        double trM0M1 = trace(M0M1).real();
                        double trM0M2 = trace(M0M2).real();
                        double trM0M3 = trace(M0M3).real();
                        double trM1M2 = trace(M1M2).real();
                        double trM1M3 = trace(M1M3).real();
                        double trM2M3 = trace(M2M3).real();
                        double trM0 = trace(mat[i[0]]).real();
                        double trM1 = trace(mat[i[1]]).real();
                        double trM2 = trace(mat[i[2]]).real();
                        double trM3 = trace(mat[i[3]]).real();
                        
                        // add to total

                        // tr4 terms
                        cx_double T1 = trM0M1M2M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM0M1M2M3);
                        res += dim*2.*(cliff*T1).real();

                        // tr3tr1 terms
                        cx_double T3 = (double)(eps[i[3]])*trM0M1M2 + (double)(eps[i[0]]*eps[i[1]]*eps[i[2]])*conj(trM0M1M2);
                        T3 = trM3*T3;

                        cx_double T4 = (double)(eps[i[2]])*trM0M1M3 + (double)(eps[i[0]]*eps[i[1]]*eps[i[3]])*conj(trM0M1M3);
                        T3 += trM2*T4;

                        cx_double T5 = (double)(eps[i[1]])*trM0M2M3 + (double)(eps[i[0]]*eps[i[2]]*eps[i[3]])*conj(trM0M2M3);
                        T3 += trM1*T5;

                        cx_double T6 = (double)(eps[i[0]])*trM1M2M3 + (double)(eps[i[1]]*eps[i[2]]*eps[i[3]])*conj(trM1M2M3);
                        T3 += trM0*T6;

                        res += 2.*(cliff*T3).real();

                        // tr2tr2 terms
                        double T7 = trM0M1*trM2M3*(eps[i[0]]*eps[i[1]] + eps[i[2]]*eps[i[3]]);
                        double T8 = trM0M2*trM1M3*(eps[i[0]]*eps[i[2]] + eps[i[1]]*eps[i[3]]);
                        double T9 = trM0M3*trM1M2*(eps[i[0]]*eps[i[3]] + eps[i[1]]*eps[i[2]]);

                        res += 2.*cliff.real()*(T7+T8+T9);
                    }
                }
            }
        }
    }

    delete [] i;
    return res;
}




double Geom24::calculate_S_old() const
{
    return g2*dirac2() + dirac4_old();
}




double Geom24::delta2(const int& x, const int& I, const int& J, const cx_double& z)
{
    if(I != J)
        return 4.*dim_omega*dim*( 2.*(z*mat[x](J,I)).real() + norm(z) );
    else
    {
        double trM = trace(mat[x]).real();
        return 8.*dim_omega*z.real()*( dim*(mat[x](I,I).real() + z.real()) + eps[x]*(trM + z.real()) );
    }
    
}

double Geom24::delta4(const int& x, const int& I, const int& J, const cx_double& z)
{
    double res = 0.;


    // D^3 dD part
    for(int i3=0; i3<nHL; ++i3)
    {
        for(int i2=0; i2<nHL; ++i2)
        {
            for(int i1=0; i1<=i3; ++i1)
            {
                cx_double cliff = omega_table_4[x + nHL*(i3 + nHL*(i2 + nHL*i1))];

                if(cliff.real() != 0. || cliff.imag() != 0.)
                {
                    // compute necessary matrix products
                    cx_mat M1M2 = mat[i1]*mat[i2];
                    cx_mat M2M3 = mat[i2]*mat[i3];
                    cx_mat M1M3 = mat[i1]*mat[i3];
                    cx_mat M1M2M3 = mat[i1]*M2M3;

                    // compute necessary traces
                    double trM1 = trace(mat[i1]).real();
                    double trM2 = trace(mat[i2]).real();
                    double trM3 = trace(mat[i3]).real();
                    double trM1M2 = trace(M1M2).real();
                    double trM2M3 = trace(M2M3).real();
                    double trM1M3 = trace(M1M3).real();
                    cx_double trM1M2M3 = trace(M1M2M3);
                    
                    // off-diagonal update
                    if(I != J)
                    {

                        // compute terms
                        // _______________________________________________________________________________________
                        cx_double T1 = M1M2M3(J,I)*z + M1M2M3(I,J)*conj(z);
                        T1 = T1 + conj(T1)*(double)(eps[i1]*eps[i2]*eps[i3]*eps[x]);
                        T1 *= (double)dim;

                        cx_double T2 = M1M2(J,I)*z + M1M2(I,J)*conj(z);
                        T2 = T2*(double)(eps[i3]) + conj(T2)*(double)(eps[i1]*eps[i2]*eps[x]);
                        T2 = T2*trM3;
                        T1 += T2;

                        cx_double T3 = M1M3(J,I)*z + M1M3(I,J)*conj(z);
                        T3 = T3*(double)(eps[i2]) + conj(T3)*(double)(eps[i1]*eps[i3]*eps[x]);
                        T3 = T3*trM2;
                        T1 += T3;

                        cx_double T4 = M2M3(J,I)*z + M2M3(I,J)*conj(z);
                        T4 = T4*(double)(eps[i1]) + conj(T4)*(double)(eps[i2]*eps[i3]*eps[x]);
                        T4 = T4*trM1;
                        T1 += T4;

                        double T5 = trM1M2*(eps[i1]*eps[i2] + eps[i3]*eps[x]);
                        T5 *= 2.*(mat[i3](J,I)*z).real();
                        T1 += T5;
                        
                        double T6 = trM2M3*(eps[i2]*eps[i3] + eps[i1]*eps[x]);
                        T6 *= 2.*(mat[i1](J,I)*z).real();
                        T1 += T6;
                        
                        double T7 = trM1M3*(eps[i1]*eps[i3] + eps[i2]*eps[x]);
                        T7 *= 2.*(mat[i2](J,I)*z).real();
                        T1 += T7;
                        //________________________________________________________________________________________
                        
                        

                        // add to total
                        if(i1 != i3)
                            res += 2.*(cliff*T1).real();
                        else
                            res += (cliff*T1).real();
                    }

                    
                    // diagonal update
                    else
                    {
                        
                        // compute terms
                        // _______________________________________________________________________________________
                        cx_double T1 = M1M2M3(I,I);
                        T1 = T1 + conj(T1)*(double)(eps[i1]*eps[i2]*eps[i3]*eps[x]);
                        T1 = T1*(double)dim;

                        cx_double T2 = M1M2(I,I);
                        T2 = T2*(double)(eps[i3]) + conj(T2)*(double)(eps[i1]*eps[i2]*eps[x]);
                        T2 *= trM3;
                        T1 += T2;

                        cx_double T3 = M1M3(I,I);
                        T3 = T3*(double)(eps[i2]) + conj(T3)*(double)(eps[i1]*eps[i3]*eps[x]);
                        T3 *= trM2;
                        T1 += T3;

                        cx_double T4 = M2M3(I,I);
                        T4 = T4*(double)(eps[i1]) + conj(T4)*(double)(eps[i2]*eps[i3]*eps[x]);
                        T4 *= trM1;
                        T1 += T4;

                        double T5 = trM1M2*(eps[i1]*eps[i2] + eps[i3]*eps[x]);
                        T5 *= mat[i3](I,I).real();
                        T1 += T5;

                        double T6 = trM2M3*(eps[i2]*eps[i3] + eps[i1]*eps[x]);
                        T6 *= mat[i1](I,I).real();
                        T1 += T6;

                        double T7 = trM1M3*(eps[i1]*eps[i3] + eps[i2]*eps[x]);
                        T7 *= mat[i2](I,I).real();
                        T1 += T7;
                        
                        cx_double T8 = conj(trM1M2M3)*(double)(eps[i1]*eps[i2]*eps[i3]) + trM1M2M3*(double)(eps[x]);
                        T1 += T8;
                        //________________________________________________________________________________________
                        
                        

                        // add to total
                        if(i1 != i3)
                            res += (cliff*T1).real()*4.*z.real();
                        else
                            res += (cliff*T1).real()*2.*z.real();
                    }
                }
            }
        }
    }

    res *= 4.;

   

    // D^2 dD^2 and D dD D dD term
    double temp = 0;
    for(int i=0; i<nHL; ++i)
    {
        double cliff = omega_table_4[x + nHL*(i + nHL*(x + nHL*i))].real();
        
        // compute necessary matrix products
        cx_mat M1M1 = mat[i]*mat[i];

        // compute necessary traces
        double trM1 = trace(mat[i]).real();
        double trM1M1 = trace(M1M1).real();

        // off-diagonal update
        if(I != J)
        {
            // compute terms D^2 dD^2
            // _______________________________________________________________________________________
            double T11 = 2*dim*( M1M1(I,I).real() + M1M1(J,J).real() );
            double T21 = 4*eps[i]*trM1*( mat[i](I,I).real() + mat[i](J,J).real() );
            double T31 = (z*mat[i](J,I)).real();
            T31 *= T31*16*eps[i]*eps[x];
            //________________________________________________________________________________________
            
            // compute terms D dD D dD
            // _______________________________________________________________________________________
            
            double T12 = (mat[i](J,I)*mat[i](J,I)*z*z).real();
            T12 += mat[i](I,I).real()*mat[i](J,J).real()*norm(z);
            T12 *= 4*dim;
            
            double T22 = 4*eps[i]*trM1*( mat[i](I,I).real() + mat[i](J,J).real() );
            double T32 = (mat[i](J,I)*z).real();
            T32 *= T32*16*eps[i]*eps[x];
            //________________________________________________________________________________________
                    
                    
            
            // add to total
            temp += 2.*dim_omega*(norm(z)*(T11+T21+4.*trM1M1) + T31);
            temp += cliff*(T12 + norm(z)*(T22+4.*trM1M1) + T32);

        } 

        // diagonal update
        else
        {
            // compute terms D^2 dD^2
            // _______________________________________________________________________________________
            double T11 = 2.*dim*M1M1(I,I).real();
            double T21 = 4.*eps[x]*M1M1(I,I).real();
            double T31 = 4.*eps[i]*trM1*mat[i](I,I).real();
            double T41 = mat[i](I,I).real();
            T41 *= T41*4.*eps[i]*eps[x];
            //________________________________________________________________________________________
                    
            // compute terms D dD D dD
            // _______________________________________________________________________________________
            double T12 = mat[i](I,I).real();
            T12 *= T12*2.*dim;
            double T22 = 4.*eps[x]*M1M1(I,I).real();
            double T32 = 4.*eps[i]*trM1*mat[i](I,I).real();
            double T42 = mat[i](I,I).real();
            T42 *= T42*4.*eps[i]*eps[x];
            //________________________________________________________________________________________
            
            // add to total
            temp += 8.*z.real()*z.real()*dim_omega*(T11+T21+T31+T41+2.*trM1M1);
            temp += 4.*z.real()*z.real()*cliff*(T12+T22+T32+T42+2.*trM1M1);
        }
    }

    res += 2.*temp;
    
                    

    // D dD^3 term
    
    // off-diangonal update
    if(I != J)
    {
        temp = 4.*dim_omega*(dim+6)*norm(z)*(mat[x](J,I)*z).real();
        res += 4.*temp;
    }

    // diagonal update
    else
    {
        double trMx = trace(mat[x]).real();
        double rez = 2.*z.real();
        temp = 2.*rez*rez*rez*dim_omega*(mat[x](I,I).real()*(dim+3.*eps[x]+3.) + eps[x]*trMx);
        res += 4.*temp;
    }


    // dD^4 term

    // off-diangonal update
    if(I != J)
    {
        temp = dim_omega*4.*norm(z)*norm(z)*(dim+6.);
        res += temp;
    }
    
    // diagonal update
    else
    {
        double rez = z.real();
        temp = dim_omega*32.*(dim+3.+4*eps[x])*rez*rez*rez*rez;
        res += temp;
    }

    return res;
}


double Geom24::delta24(const int& x, const int& I, const int& J, const cx_double& z)
{
    return g2*delta2(x,I,J,z)+delta4(x,I,J,z);
}



void Geom24::delta24_debug(const double& scale, const int& iter, gsl_rng* engine, ostream& out_s)
{
    // iter iterations
    for(int i=0; i<iter; ++i)
    {
        // find Si
        shuffle(engine);
        //double Si_1 = dirac2();
        //double Si_1 = dirac4();
        double Si_1 = calculate_S();
        
        // aandom move
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

        //double dS = delta2(x, I, J, z);
        //double dS = delta4(x, I, J, z);
        double dS = delta24(x, I, J, z);

        bool diag = true;
        // update matrix element
        if(I != J)
        {
            mat[x](I,J) += z;
            mat[x](J,I) += conj(z);
            diag = false;
        }
        else
            mat[x](I,I) += 2.*z;
        
        // find Sf
        //double Sf_1 = dirac2();
        //double Sf_1 = dirac4();
        double Sf_1 = calculate_S();

    
        // print S
        out_s.precision(16);
        out_s << Si_1 << " " << Sf_1 << " " << Sf_1-Si_1 << " " << dS << " " << Sf_1-Si_1-dS << " " << diag << endl; 
    }
}
