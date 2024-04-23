// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues and Manuel Assunção.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @metric_covariant.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_COVARIANT
#define GYRONIMO_METRIC_COVARIANT

#include <gyronimo/core/contraction.hh>
#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/core/SM3algebra.hh>

namespace gyronimo {

//! Abstract class implementing basic 3x3 metric-tensor functionality.
/*!
    Requires derived classes to implement the 6 independent components @f$
    g_{ij} @f$ of the covariant metric tensor for a given coordinate system in
    @f$\mathbb{R}^3@f$ [via the `operator()()` method] and their corresponding
    18 partial derivatives [via the `del()` method]. The class implements
    general methods to compute the Jacobian (square root of its determinant) and
    its gradient, to transform contravariant vectors into covariant ones and
    conversely, to get its inverse derivatives, and to get Christoffel's
    symbols. These methods are left as virtual to allow more efficient
    reimplementations in derived classes, if needed. The `IR3 q` argument stands
    for the three independent variables mapping the underlying coordinate system
    to the cartesian space via some invertible mapping @f$ \mathbf{x}(q^u, q^v,
    q^w) @f$. Notice that such mapping (i.e., a `gyronimo::morphism`) does not
    have to be explicitly known at this level, only the components @f$ g_{ij} =
    \partial_i \mathbf{x} \cdot \partial_j \mathbf{x} @f$ and their derivatives.
    Regarding units, all methods must ensure that @f$g_{ij} q^i q^j@f$ returns
    values in SI (m^2). **Atention**: the coordinates @f$ \{q^u, q^v, q^w\} @f$
    must be right-hand ordered to ensure a positive determinant.
*/
class metric_covariant {
 public:
  metric_covariant() {};
  virtual ~metric_covariant() {};

  virtual SM3 operator()(const IR3& q) const = 0;
  virtual dSM3 del(const IR3& q) const = 0;

  virtual double jacobian(const IR3& q) const;
  virtual IR3 del_jacobian(const IR3& q) const;
  virtual IR3 to_covariant(const IR3& B, const IR3& q) const;
  virtual IR3 to_contravariant(const IR3& B, const IR3& q) const;
  virtual SM3 inverse(const IR3& q) const;
  virtual dSM3 del_inverse(const IR3& q) const;
  virtual ddIR3 christoffel_first_kind(const IR3& q) const;
  virtual ddIR3 christoffel_second_kind(const IR3& q) const;
  virtual IR3 inertial_force(const IR3& q, const IR3& dot_q) const;
};

inline SM3 metric_covariant::inverse(const IR3& q) const {
  SM3 m = (*this)(q);
  return gyronimo::inverse(m);
}

//! General index-lowering of contravariant vector `B` at position `q`.
inline IR3 metric_covariant::to_covariant(const IR3& B, const IR3& q) const {
  return contraction((*this)(q), B);
}

//! General index-raising of covariant vector `B` at position `q`.
inline IR3 metric_covariant::to_contravariant(
    const IR3& B, const IR3& q) const {
  return contraction(this->inverse(q), B);
}

//! General Christoffel symbol @f$\Gamma^k_{ij} = g^{km} \, \Gamma_{mij}@f$.
inline ddIR3 metric_covariant::christoffel_second_kind(const IR3& q) const {
  return contraction<first>(
      this->christoffel_first_kind(q), this->inverse(q));
}

}  // namespace gyronimo

#endif  // GYRONIMO_METRIC_COVARIANT
