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

// @morphism.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM
#define GYRONIMO_MORPHISM

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/core/contraction.hh>

namespace gyronimo {

//! Abstract morphism from curvilinear `q` into *cartesian* `x` coordinates.
/*!
    Derived classes must implement the map @f$ \mathbf{x}(q^\alpha) @f$ from the
    curvilinear coordinates @f$ q^\alpha @f$ into the cartesian space and its
    inverse, along with the first-order and second-order partial derivatives of
    the former. The default implementation provides the inverse derivative and
    the transformation jacobian, dual and tangent-space basis, and conversion
    between covariant and contravariant components. These methods are left
    virtual to allow more efficient reimplementations in derived classes, if
    needed.
*/
class morphism {
 public:
  morphism() {};
  virtual ~morphism() {};

  virtual IR3 operator()(const IR3& q) const = 0;
  virtual IR3 inverse(const IR3& x) const = 0;
  virtual dIR3 del(const IR3& q) const = 0;
  virtual ddIR3 ddel(const IR3& q) const = 0;

  virtual double jacobian(const IR3& q) const;
  virtual dIR3 del_inverse(const IR3& q) const;
  virtual dIR3 tan_basis(const IR3& q) const;
  virtual dIR3 dual_basis(const IR3& q) const;
  virtual IR3 to_covariant(const IR3& A, const IR3& q) const;
  virtual IR3 to_contravariant(const IR3& A, const IR3& q) const;
  virtual IR3 from_covariant(const IR3& A, const IR3& q) const;
  virtual IR3 from_contravariant(const IR3& A, const IR3& q) const;
  virtual IR3 translation(const IR3& q, const IR3& delta) const;
};

}  // end namespace gyronimo

#endif  // GYRONIMO_MORPHISM
