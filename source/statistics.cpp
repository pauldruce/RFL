#include <armadillo>
#include <cmath>
#include "statistics.hpp"

using namespace std;
using namespace arma;


void jackknife(const vec& vec_uncorr, double& avg, double& var, double f(const vec&))
{
    // Find vector size
    int size = vec_uncorr.n_elem;

    // Create vector of delete-1 clusters
    vec vec_del1(size);
    
    // Calculate delete-1 clusters
    for(int i=0; i<size; ++i)
    {
        // Create i-th cluster
        vec ith_cluster(size-1);
        for(int j=0; j<size; ++j)
        {
            // Copy element in the same position if j < i
            if(j < i)
                ith_cluster(j) = vec_uncorr(j);
            // Copy element one position behind if j > i
            else if(j > i)
                ith_cluster(j-1) = vec_uncorr(j);

        }

        // Put the return value of f into vec_del1 in i-th position
        vec_del1(i) = f(ith_cluster);
    }

    // Calculate mean
    avg = mean(vec_del1);

    // Calculate variance
    var = 0;
    for(int i=0; i<size; ++i)
        var += pow(vec_del1(i) - avg, 2);
    var *= (double)(size-1)/size;
}

double my_mean(const vec& vec_uncorr)
{
    return mean(vec_uncorr);
}

double my_var(const vec& vec_uncorr)
{
    return var(vec_uncorr, 1);
}

double my_sus(const vec& vec_uncorr)
{
    return var(vec_uncorr, 1)/mean(vec_uncorr);
}
