//
// Created by Paul Druce on 15/04/2022.
//
#include "doctest/doctest.h"
#include "gsl/gsl_rng.h"
#include "geometry.hpp"
#include <armadillo>

using namespace arma;

static auto engine = gsl_rng_alloc(gsl_rng_ranlxd1);
constexpr int N = 10;

static void CompareDeltaToDifference(int p, int q, int dim, double g2)
{
    Geom24 G(p, q, dim, g2);

    for (int i = 0; i < N; ++i)
    {
        G.shuffle(engine);
        double Si = G.calculate_S();

        // random entry
        int x = (int)(G.get_nHL() * gsl_rng_uniform(engine));
        int I = (int)(G.get_dim() * gsl_rng_uniform(engine));
        int J = (int)(G.get_dim() * gsl_rng_uniform(engine));
        double re;
        double im;
        cx_double z;
        if (I != J)
        {
            re = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
            im = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
            z = cx_double(re, im);
        }
        else
        {
            re = 0.5 * (-1. + 2. * gsl_rng_uniform(engine));
            z = cx_double(re, 0);
        }

        double dS = G.delta24(x, I, J, z);

        cx_mat m = G.get_mat(x);
        m(I, J) += z;
        m(J, I) += conj(z);
        G.set_mat(x, m);

        double Sf = G.calculate_S();
        CHECK_MESSAGE(fabs((Sf - Si) - dS) < 1e-8,
                      "     Sf-Si =", Sf - Si, ", dS = ", dS, " and |(Sf-Si) - dS| = ", fabs(Sf - Si - dS));
    }
}

TEST_SUITE("Action Calculation Tests")
{
    TEST_CASE("Calculate the dS/Delta tests")
    {
        gsl_rng_set(engine, time(nullptr));

        // Bad way to implement parametrised tests.
        for (int dim = 3; dim < 6; dim++)
        {
            for (int p = 1; p < 4; p++)
            {
                for (int q = 1; q < 4; q++)
                {
                    CAPTURE(p);
                    CAPTURE(q);
                    CAPTURE(dim);
                    CompareDeltaToDifference(p, q, dim, -2.2);
                }
            }
        }
    }
}