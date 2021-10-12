// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @SM3algebra.hh

#ifndef GYRONIMO_SM3ALGEBRA
#define GYRONIMO_SM3ALGEBRA

#include <array>

namespace gyronimo {

//! Symmetric 3x3 matrix (i.e., SM3).
/*!
    *Aggregate* type, natively supporting list initialization. The index
    operator `B[SM3::ij]` returns the component @f$B^{ij}@f$ (for contravariant
    tensors) or @f$B_{ij}@f$ (for covariant tensors) with `i,j = u, v, w`.
*/
class SM3 {
 public:
  enum index : size_t {uu = 0, uv = 1, uw = 2, vv = 3, vw = 4, ww = 5};

  double& operator[](index i) {return data_[i];};
  const double& operator[](index i) const {return data_[i];};

  std::array<double, 6> data_;
};

//! Partial derivatives of a symmetric 3x3 matrix.
/*!
    *Aggregate* type, natively supporting list initialization. The index
    operator `B[dSM3::ijk]` returns the value @f$\partial_k B^{ij}@f$ (for
    contravariant tensors) or @f$\partial_k B_{ij}@f$ (for covariant tensors)
    with `i,j,k = u, v, w`.
*/
class dSM3 {
 public:
  enum index : size_t {
    uuu = 0, uuv = 1, uuw = 2, uvu = 3, uvv = 4, uvw = 5,
    uwu = 6, uwv = 7, uww = 8, vvu = 9, vvv = 10, vvw = 11,
    vwu = 12, vwv = 13, vww = 14, wwu = 15, wwv = 16, www = 17};

  double& operator[](index i) {return data_[i];};
  const double& operator[](index i) const {return data_[i];};

  std::array<double, 18> data_;
};

} // end namespace gyronimo.

#endif // end GYRONIMO_SM3ALGEBRA
