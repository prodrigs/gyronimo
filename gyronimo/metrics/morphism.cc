// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Paulo Rodrigues and Manuel Assunção.

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

// @morphism.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

dIR3 morphism::del_inverse(const IR3& q) const {
  return gyronimo::inverse(del(q));
}

//! Tangent basis @f$ \textbf{e}_\alpha = \partial_\alpha \mathbf{x} @f$.
dIR3 morphism::tan_basis(const IR3& q) const { return del(q); }

//! Dual basis @f$ \textbf{e}^\alpha = \nabla q^\alpha(\mathbf{x}) @f$.
dIR3 morphism::dual_basis(const IR3& q) const { return del_inverse(q); }

//! Covariant components of *cartesian* `A` at curvilinear @f$q^\alpha@f$.
IR3 morphism::to_covariant(const IR3& A, const IR3& q) const {
  return contraction<first>(del(q), A);
}

//! Contravariant components of *cartesian* `A` at curvilinear @f$q^\alpha@f$.
IR3 morphism::to_contravariant(const IR3& A, const IR3& q) const {
  return contraction<second>(del_inverse(q), A);
}

//! Cartesian vector from covariant components.
IR3 morphism::from_covariant(const IR3& A, const IR3& q) const {
  return contraction<first>(del_inverse(q), A);
}

//! Cartesian vector from contravariant components.
IR3 morphism::from_contravariant(const IR3& A, const IR3& q) const {
  return contraction<second>(del(q), A);
}

//! Jacobian @f$ \mathbf{e}_u \cdot (\mathbf{e}_v \times \mathbf{e}_w) @f$.
double morphism::jacobian(const IR3& q) const {
  dIR3 e = del(q);
  return e[dIR3::uu] * (e[dIR3::vv] * e[dIR3::ww] - e[dIR3::vw] * e[dIR3::wv]) +
      e[dIR3::uv] * (e[dIR3::vw] * e[dIR3::wu] - e[dIR3::vu] * e[dIR3::ww]) +
      e[dIR3::uw] * (e[dIR3::vu] * e[dIR3::wv] - e[dIR3::vv] * e[dIR3::wu]);
}

//! Curvilinear position of *cartesian* @f$ \mathbf{x}(q^\alpha) + \Delta @f$.
IR3 morphism::translation(const IR3& q, const IR3& delta) const {
  return this->inverse((*this)(q) + delta);
}

}  // end namespace gyronimo
