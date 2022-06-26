//
// Created by Paul Druce on 15/04/2022.
//
#include "doctest/doctest.h"
#include "geometry.hpp"


TEST_SUITE_BEGIN("Compilation Tests");

TEST_CASE("")
{
   int p, q, dim;
   double g2;
   dim = 5;
   for (p = 1; p < 3; p++)
   {
      for (q = 1; q < 3; q++)
      {
         for (g2 = -5.0; g2 < 0.0; g2 += 1.0)
         {
            CAPTURE(p);
            CAPTURE(q);
            CAPTURE(dim);
            CAPTURE(g2);
            CHECK_NOTHROW_MESSAGE(Geom24 G(p, q, dim, g2),
                                  "Constructing Geom24 class with p=", p, " q=", q, " dim=", dim, " and g2=", g2,
                                  "\n");
         }
      }
   }
}

TEST_SUITE_END;