//
// Created by Paul Druce on 14/04/2022.
//

#include "doctest/doctest.h"
#include <gsl/gsl_rng.h>
#include "geometry.hpp"

static auto engine = gsl_rng_alloc(gsl_rng_ranlxd1);

static void CompareActions(int p, int q, int dim, double g2)
{
   constexpr int N = 1000;
   Geom24 G(p, q, dim, g2);
   for (int i = 0; i < N; ++i)
   {
      G.shuffle(engine);
      double d = G.get_dim();

      double S1 = G.calculate_S() / (d * d);
      double S2 = G.calculate_S_from_dirac() / (d * d);
      CAPTURE(p);
      CAPTURE(q);
      CAPTURE(dim);
      CAPTURE(g2);
      CHECK_MESSAGE(fabs(S1 - S2) < 1e-8, "Methods differ more then 1e-8");
   }
}


TEST_SUITE("Action Calculation Tests")
{
   TEST_CASE("Two methods produce equal results")
   {
      gsl_rng_set(engine, time(nullptr));

      SUBCASE("Type (1,1), dim = 5, g2 = -3.0")
      {
         CompareActions(1, 1, 5, -3.0);
      }

      SUBCASE("Type (2,2), dim = 6, g2 = -2.2")
      {
         CompareActions(2, 2, 6, -2.2);
      }

      SUBCASE("Type (1,2), dim = 7, g2 = -0.5")
      {
         CompareActions(1, 2, 7, -0.5);
      }

      SUBCASE("Type (2,1), dim  4, g2 = -0.1")
      {
         CompareActions(2, 1, 4, -0.1);
      }

      SUBCASE("Type (0,5), dim = 6, g2 = -2.8")
      {
         CompareActions(0, 5, 4, -2.8);
      }
   }
}


