#include "Geom24.hpp"

using namespace std;
using namespace arma;

// PD: All methods moved into the DiracOperator class.

cx_mat Geom24::der_dirac24(const int& k, const bool& herm) const {
  return g2 * der_dirac2(k) + der_dirac4(k, herm);
}

cx_mat Geom24::der_dirac2(const int& k) const {
  cx_mat res(dim, dim, fill::eye);

  res *= eps[k] * trace(mat[k]).real();
  res += dim * mat[k];

  return 4 * dim_omega * res;
}

cx_mat Geom24::der_dirac4(const int& k, const bool& herm) const {
  cx_mat res(dim, dim, fill::zeros);

  // four distinct indices
  for (int i1 = 0; i1 < nHL; ++i1) {
    if (i1 != k) {
      for (int i2 = i1 + 1; i2 < nHL; ++i2) {
        if (i2 != k) {
          for (int i3 = i2 + 1; i3 < nHL; ++i3) {
            if (i3 != k) {
              // epsilon factor
              double e = eps[k] * eps[i1] * eps[i2] * eps[i3];

              if (e < 0) {
                // clifford product
                double cliff1 = omega_table_4[i3 + nHL * (i2 + nHL * (i1 + nHL * k))].imag();
                double cliff2 = omega_table_4[i2 + nHL * (i3 + nHL * (i1 + nHL * k))].imag();
                double cliff3 = omega_table_4[i3 + nHL * (i1 + nHL * (i2 + nHL * k))].imag();

                if (fabs(cliff1) > 1e-10) {
                  res += compute_B4(k, i1, i2, i3, cliff1, true);
                  res += compute_B4(k, i1, i3, i2, cliff2, true);
                  res += compute_B4(k, i2, i1, i3, cliff3, true);
                }
              } else {
                // clifford product
                double cliff1 = omega_table_4[i3 + nHL * (i2 + nHL * (i1 + nHL * k))].real();
                double cliff2 = omega_table_4[i2 + nHL * (i3 + nHL * (i1 + nHL * k))].real();
                double cliff3 = omega_table_4[i3 + nHL * (i1 + nHL * (i2 + nHL * k))].real();

                if (fabs(cliff1) > 1e-10) {
                  res += compute_B4(k, i1, i2, i3, cliff1, false);
                  res += compute_B4(k, i1, i3, i2, cliff2, false);
                  res += compute_B4(k, i2, i1, i3, cliff3, false);
                }
              }
            }
          }
        }
      }
    }
  }
  res = res + res.t();

  // two distinct pairs of equal indices
  for (int i = 0; i < nHL; ++i) {
    if (i != k) {
      res += compute_B2(k, i);
    }
  }

  // all indices equal
  res += compute_B(k);

  if (herm) {
    return 2 * (res + res.t());
  } else {
    return 4 * res;
  }
}

cx_mat Geom24::compute_B4(const int& k,
                          const int& i2,
                          const int& i3,
                          const int& i4,
                          const double& cliff,
                          const bool& neg) const {
  // base matrix products
  cx_mat M2M3 = mat[i2] * mat[i3];
  cx_mat M2M4 = mat[i2] * mat[i4];
  cx_mat M3M4 = mat[i3] * mat[i4];
  cx_mat M2M3M4 = M2M3 * mat[i4];

  // return value
  cx_mat res(dim, dim, fill::eye);

  if (neg) {

    // traces
    double tr234 = trace(M2M3M4).imag();
    double tr2 = trace(mat[i2]).real();
    double tr3 = trace(mat[i3]).real();
    double tr4 = trace(mat[i4]).real();

    // compute sum
    res *= -2 * eps[k] * tr234;
    res += cx_double(0., dim) * (M2M3M4 - M2M3M4.t());
    res += cx_double(0., eps[i2] * tr2) * (M3M4 - M3M4.t());
    res += cx_double(0., eps[i3] * tr3) * (M2M4 - M2M4.t());
    res += cx_double(0., eps[i4] * tr4) * (M2M3 - M2M3.t());
  } else {
    // traces
    double tr234 = trace(M2M3M4).real();
    double tr23 = trace(M2M3).real();
    double tr24 = trace(M2M4).real();
    double tr34 = trace(M3M4).real();
    double tr2 = trace(mat[i2]).real();
    double tr3 = trace(mat[i3]).real();
    double tr4 = trace(mat[i4]).real();

    // compute sum
    res *= 2 * eps[k] * tr234;
    res += dim * (M2M3M4 + M2M3M4.t());
    res += eps[i2] * tr2 * (M3M4 + M3M4.t());
    res += eps[i3] * tr3 * (M2M4 + M2M4.t());
    res += eps[i4] * tr4 * (M2M3 + M2M3.t());
    res += 2 * eps[k] * eps[i2] * tr34 * mat[i2];
    res += 2 * eps[k] * eps[i3] * tr24 * mat[i3];
    res += 2 * eps[k] * eps[i4] * tr23 * mat[i4];
  }

  return cliff * res;
}

cx_mat Geom24::compute_B2(const int& k, const int& i) const {
  // clifford product
  double cliff = omega_table_4[i + nHL * (k + nHL * (i + nHL * k))].real();

  // base matrix products
  cx_mat MiMk = mat[i] * mat[k];
  cx_mat MiMi = mat[i] * mat[i];
  cx_mat MiMiMk = mat[i] * MiMk;
  cx_mat MiMkMi = MiMk * mat[i];

  // traces
  double triki = trace(MiMkMi).real();
  double trik = trace(MiMk).real();
  double trii = trace(MiMi).real();
  double tri = trace(mat[i]).real();
  double trk = trace(mat[k]).real();

  // return value
  cx_mat res(dim, dim, fill::eye);

  if (cliff < 0) {
    // compute sum
    res *= eps[k] * triki;
    res += dim * (MiMiMk + MiMiMk.t() - MiMkMi);
    res += eps[i] * tri * (MiMk + MiMk.t());
    res += 2 * eps[k] * eps[i] * trik * mat[i];
    res += eps[k] * trk * MiMi;
    res += trii * mat[k];
  } else {
    // compute sum
    res *= 3 * eps[k] * triki;
    res += dim * (MiMiMk + MiMiMk.t() + MiMkMi);
    res += 3 * eps[i] * tri * (MiMk + MiMk.t());
    res += 6 * eps[k] * eps[i] * trik * mat[i];
    res += 3 * eps[k] * trk * MiMi;
    res += 3 * trii * mat[k];
  }

  return 2 * dim_omega * res;
}

cx_mat Geom24::compute_B(const int& k) const {
  // base matrix products
  cx_mat M2 = mat[k] * mat[k];
  cx_mat M3 = mat[k] * M2;

  // traces
  double tr3 = trace(M3).real();
  double tr2 = trace(M2).real();
  double tr1 = trace(mat[k]).real();

  cx_mat res(dim, dim, fill::eye);
  res *= eps[k] * tr3;
  res += dim * M3;
  res += 3 * tr2 * mat[k];
  res += 3 * eps[k] * tr1 * M2;

  return 2 * dim_omega * res;
}
