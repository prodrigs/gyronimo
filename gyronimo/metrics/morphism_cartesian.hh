// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Manuel Assunção and Paulo Rodrigues.

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

// @morphism_cartesian.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_CARTESIAN
#define GYRONIMO_MORPHISM_CARTESIAN

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Trivial (identity) morphism from cartesian coordinates (in SI units).
class morphism_cartesian : public morphism {
 public:
  morphism_cartesian() : morphism() {};
  virtual ~morphism_cartesian() override {};

  virtual IR3 operator()(const IR3& q) const override { return q; };
  virtual IR3 inverse(const IR3& x) const override { return x; };
  virtual dIR3 del(const IR3& q) const override;
  virtual ddIR3 ddel(const IR3& q) const override;

  virtual double jacobian(const IR3& q) const override { return 1; };
  virtual dIR3 del_inverse(const IR3& q) const override;
  virtual IR3 to_covariant(const IR3& A, const IR3& q) const override;
  virtual IR3 to_contravariant(const IR3& A, const IR3& q) const override;
  virtual IR3 from_covariant(const IR3& A, const IR3& q) const override;
  virtual IR3 from_contravariant(
      const IR3& A, const IR3& q) const override;
  virtual IR3 translation(const IR3& q, const IR3& delta) const override;
};

inline IR3 morphism_cartesian::translation(
    const IR3& q, const IR3& delta) const {
  return q + delta;
}
inline dIR3 morphism_cartesian::del(const IR3& q) const {
  return {1, 0, 0, 0, 1, 0, 0, 0, 1};
}
inline ddIR3 morphism_cartesian::ddel(const IR3& q) const {
  return {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}
inline dIR3 morphism_cartesian::del_inverse(const IR3& q) const {
  return {1, 0, 0, 0, 1, 0, 0, 0, 1};
}
inline IR3 morphism_cartesian::to_covariant(const IR3& A, const IR3& q) const {
  return A;
}
inline IR3 morphism_cartesian::to_contravariant(
    const IR3& A, const IR3& q) const {
  return A;
}
inline IR3 morphism_cartesian::from_covariant(
    const IR3& A, const IR3& q) const {
  return A;
}
inline IR3 morphism_cartesian::from_contravariant(
    const IR3& A, const IR3& q) const {
  return A;
}

}  // end namespace gyronimo

#endif  // GYRONIMO_MORPHISM_CARTESIAN
