#ifndef STATS_HPP
#define STATS_HPP

#include <armadillo>

// Jack knife mean and variance estimate of an arbitrary function f
void jackknife(const arma::vec&, double&, double&, double f(const arma::vec&));

// Some stat functions
double my_mean(const arma::vec&);
double my_var(const arma::vec&);
double my_sus(const arma::vec&);

#endif

