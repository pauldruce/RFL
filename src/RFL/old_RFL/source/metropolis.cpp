#include "Geom24.hpp"
#include <cmath>

using namespace std;
using namespace arma;

// PD: Moved into Metropolis class for now.

double Geom24::delta24(const int& x, const int& I, const int& J, const cx_double& z) {
  return g2 * delta2(x, I, J, z) + delta4(x, I, J, z);
}

double Geom24::delta2(const int& x, const int& I, const int& J, const cx_double& z) {
  if (I != J) {
    return 4. * dim_omega * dim * (2. * (z * mat[x](J, I)).real() + norm(z));
  } else {
    double trM = trace(mat[x]).real();
    return 8. * dim_omega * z.real() * (dim * (mat[x](I, I).real() + z.real()) + eps[x] * (trM + z.real()));
  }
}

double Geom24::delta4(const int& x, const int& I, const int& J, const cx_double& z) {
  double res = 0.;

  // D^3 dD part
  for (int i3 = 0; i3 < nHL; ++i3) {
    for (int i2 = 0; i2 < nHL; ++i2) {
      for (int i1 = 0; i1 <= i3; ++i1) {
        cx_double cliff = omega_table_4[x + nHL * (i3 + nHL * (i2 + nHL * i1))];

        if (fabs(cliff.real()) > 1e-10 || fabs(cliff.imag()) > 1e-10) {
          // compute necessary matrix products
          cx_mat M1M2 = mat[i1] * mat[i2];
          cx_mat M2M3 = mat[i2] * mat[i3];
          cx_mat M1M3 = mat[i1] * mat[i3];
          cx_mat M1M2M3 = mat[i1] * M2M3;

          // compute necessary traces
          double trM1 = trace(mat[i1]).real();
          double trM2 = trace(mat[i2]).real();
          double trM3 = trace(mat[i3]).real();
          double trM1M2 = trace(M1M2).real();
          double trM2M3 = trace(M2M3).real();
          double trM1M3 = trace(M1M3).real();
          cx_double trM1M2M3 = trace(M1M2M3);

          // off-diagonal update
          if (I != J) {

            // compute terms
            // _______________________________________________________________________________________
            cx_double T1 = M1M2M3(J, I) * z + M1M2M3(I, J) * conj(z);
            T1 = T1 + conj(T1) * (double)(eps[i1] * eps[i2] * eps[i3] * eps[x]);
            T1 *= (double)dim;

            cx_double T2 = M1M2(J, I) * z + M1M2(I, J) * conj(z);
            T2 = T2 * (double)(eps[i3]) + conj(T2) * (double)(eps[i1] * eps[i2] * eps[x]);
            T2 = T2 * trM3;
            T1 += T2;

            cx_double T3 = M1M3(J, I) * z + M1M3(I, J) * conj(z);
            T3 = T3 * (double)(eps[i2]) + conj(T3) * (double)(eps[i1] * eps[i3] * eps[x]);
            T3 = T3 * trM2;
            T1 += T3;

            cx_double T4 = M2M3(J, I) * z + M2M3(I, J) * conj(z);
            T4 = T4 * (double)(eps[i1]) + conj(T4) * (double)(eps[i2] * eps[i3] * eps[x]);
            T4 = T4 * trM1;
            T1 += T4;

            double T5 = trM1M2 * (eps[i1] * eps[i2] + eps[i3] * eps[x]);
            T5 *= 2. * (mat[i3](J, I) * z).real();
            T1 += T5;

            double T6 = trM2M3 * (eps[i2] * eps[i3] + eps[i1] * eps[x]);
            T6 *= 2. * (mat[i1](J, I) * z).real();
            T1 += T6;

            double T7 = trM1M3 * (eps[i1] * eps[i3] + eps[i2] * eps[x]);
            T7 *= 2. * (mat[i2](J, I) * z).real();
            T1 += T7;
            //________________________________________________________________________________________

            // add to total
            if (i1 != i3) {
              res += 2. * (cliff * T1).real();
            } else {
              res += (cliff * T1).real();
            }
          }

          // diagonal update
          else {

            // compute terms
            // _______________________________________________________________________________________
            cx_double T1 = M1M2M3(I, I);
            T1 = T1 + conj(T1) * (double)(eps[i1] * eps[i2] * eps[i3] * eps[x]);
            T1 = T1 * (double)dim;

            cx_double T2 = M1M2(I, I);
            T2 = T2 * (double)(eps[i3]) + conj(T2) * (double)(eps[i1] * eps[i2] * eps[x]);
            T2 *= trM3;
            T1 += T2;

            cx_double T3 = M1M3(I, I);
            T3 = T3 * (double)(eps[i2]) + conj(T3) * (double)(eps[i1] * eps[i3] * eps[x]);
            T3 *= trM2;
            T1 += T3;

            cx_double T4 = M2M3(I, I);
            T4 = T4 * (double)(eps[i1]) + conj(T4) * (double)(eps[i2] * eps[i3] * eps[x]);
            T4 *= trM1;
            T1 += T4;

            double T5 = trM1M2 * (eps[i1] * eps[i2] + eps[i3] * eps[x]);
            T5 *= mat[i3](I, I).real();
            T1 += T5;

            double T6 = trM2M3 * (eps[i2] * eps[i3] + eps[i1] * eps[x]);
            T6 *= mat[i1](I, I).real();
            T1 += T6;

            double T7 = trM1M3 * (eps[i1] * eps[i3] + eps[i2] * eps[x]);
            T7 *= mat[i2](I, I).real();
            T1 += T7;

            cx_double T8 =
                conj(trM1M2M3) * (double)(eps[i1] * eps[i2] * eps[i3]) + trM1M2M3 * (double)(eps[x]);
            T1 += T8;
            //________________________________________________________________________________________

            // add to total
            if (i1 != i3) {
              res += (cliff * T1).real() * 4. * z.real();
            } else {
              res += (cliff * T1).real() * 2. * z.real();
            }
          }
        }
      }
    }
  }

  res *= 4.;

  // D^2 dD^2 and D dD D dD term
  double temp = 0;
  for (int i = 0; i < nHL; ++i) {
    double cliff = omega_table_4[x + nHL * (i + nHL * (x + nHL * i))].real();

    // compute necessary matrix products
    cx_mat M1M1 = mat[i] * mat[i];

    // compute necessary traces
    double trM1 = trace(mat[i]).real();
    double trM1M1 = trace(M1M1).real();

    // off-diagonal update
    if (I != J) {
      // compute terms D^2 dD^2
      // _______________________________________________________________________________________
      double T11 = 2 * dim * (M1M1(I, I).real() + M1M1(J, J).real());
      double T21 = 4 * eps[i] * trM1 * (mat[i](I, I).real() + mat[i](J, J).real());
      double T31 = (z * mat[i](J, I)).real();
      T31 *= T31 * 16 * eps[i] * eps[x];
      //________________________________________________________________________________________

      // compute terms D dD D dD
      // _______________________________________________________________________________________

      double T12 = (mat[i](J, I) * mat[i](J, I) * z * z).real();
      T12 += mat[i](I, I).real() * mat[i](J, J).real() * norm(z);
      T12 *= 4 * dim;

      double T22 = 4 * eps[i] * trM1 * (mat[i](I, I).real() + mat[i](J, J).real());
      double T32 = (mat[i](J, I) * z).real();
      T32 *= T32 * 16 * eps[i] * eps[x];
      //________________________________________________________________________________________

      // add to total
      temp += 2. * dim_omega * (norm(z) * (T11 + T21 + 4. * trM1M1) + T31);
      temp += cliff * (T12 + norm(z) * (T22 + 4. * trM1M1) + T32);
    }

    // diagonal update
    else {
      // compute terms D^2 dD^2
      // _______________________________________________________________________________________
      double T11 = 2. * dim * M1M1(I, I).real();
      double T21 = 4. * eps[x] * M1M1(I, I).real();
      double T31 = 4. * eps[i] * trM1 * mat[i](I, I).real();
      double T41 = mat[i](I, I).real();
      T41 *= T41 * 4. * eps[i] * eps[x];
      //________________________________________________________________________________________

      // compute terms D dD D dD
      // _______________________________________________________________________________________
      double T12 = mat[i](I, I).real();
      T12 *= T12 * 2. * dim;
      double T22 = 4. * eps[x] * M1M1(I, I).real();
      double T32 = 4. * eps[i] * trM1 * mat[i](I, I).real();
      double T42 = mat[i](I, I).real();
      T42 *= T42 * 4. * eps[i] * eps[x];
      //________________________________________________________________________________________

      // add to total
      temp += 8. * z.real() * z.real() * dim_omega * (T11 + T21 + T31 + T41 + 2. * trM1M1);
      temp += 4. * z.real() * z.real() * cliff * (T12 + T22 + T32 + T42 + 2. * trM1M1);
    }
  }

  res += 2. * temp;

  // D dD^3 term

  // off-diagonal update
  if (I != J) {
    temp = 4. * dim_omega * (dim + 6) * norm(z) * (mat[x](J, I) * z).real();
    res += 4. * temp;
  }

  // diagonal update
  else {
    double trMx = trace(mat[x]).real();
    double rez = 2. * z.real();
    temp = 2. * rez * rez * rez * dim_omega * (mat[x](I, I).real() * (dim + 3. * eps[x] + 3.) + eps[x] * trMx);
    res += 4. * temp;
  }

  // dD^4 term

  // off-diagonal update
  if (I != J) {
    temp = dim_omega * 4. * norm(z) * norm(z) * (dim + 6.);
    res += temp;
  }

  // diagonal update
  else {
    double rez = z.real();
    temp = dim_omega * 32. * (dim + 3. + 4 * eps[x]) * rez * rez * rez * rez;
    res += temp;
  }

  return res;
}
