//
// Created by Paul Druce on 14/04/2022.
//

#include "doctest.h"
#include <gsl/gsl_rng.h>
#include "geometry.hpp"

TEST_SUITE("Test")
{
   TEST_CASE("Basics")
   {
      gsl_rng *engine = gsl_rng_alloc(gsl_rng_ranlxd1);
      gsl_rng_set(engine, time(nullptr));


      // create geometry
      int p = 1;
      int q = 1;
      int dim = 5;
      double g2 = -3.0;
      Geom24 G(p, q, dim, g2);
      constexpr int N = 1000;
      for (int i = 0; i < N; ++i)
      {
         G.shuffle(engine);
         double d = G.get_dim();

         double S1 = G.calculate_S() / (d * d);
         double S2 = G.calculate_S_from_dirac() / (d * d);

         std::cout.precision(16);
         CHECK(fabs(S1 - S2) < 1e-8);
         if (fabs(S1 - S2) > 1e-8)
         {
            std::cout << S1 << " " << S2 << " " << S1 - S2 << std::endl;
         }
      }
   }
}


