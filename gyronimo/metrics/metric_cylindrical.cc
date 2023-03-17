// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Manuel Assunção.

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

// @metric_cylindrical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_cylindrical.hh>

#include <cmath>

namespace gyronimo {

metric_cylindrical::metric_cylindrical(const morphism_cylindrical* morph)
    : metric_connected(morph), Lref_(morph->Lref()), Lref_2_(Lref_ * Lref_),
      iLref_2_(1 / Lref_2_), Lref_3_(Lref_ * Lref_ * Lref_) {}
SM3 metric_cylindrical::operator()(const IR3& q) const {
  return {Lref_2_, 0.0, 0.0, Lref_2_ * q[IR3::u] * q[IR3::u], 0.0, Lref_2_};
}
SM3 metric_cylindrical::inverse(const IR3& q) const {
  return {
      iLref_2_, 0.0, 0.0, iLref_2_ / (q[IR3::u] * q[IR3::u]), 0.0, iLref_2_};
}
dSM3 metric_cylindrical::del(const IR3& q) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      2 * q[IR3::u] * Lref_2_, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
double metric_cylindrical::jacobian(const IR3& q) const {
  return Lref_3_ * q[IR3::u];
}
IR3 metric_cylindrical::del_jacobian(const IR3& q) const {
  return {Lref_3_, 0.0, 0.0};
}
IR3 metric_cylindrical::to_covariant(const IR3& B, const IR3& q) const {
  return {Lref_2_ * B[IR3::u],
      Lref_2_ * B[IR3::v] * q[IR3::u] * q[IR3::u], Lref_2_ * B[IR3::w]};
}
IR3 metric_cylindrical::to_contravariant(const IR3& B, const IR3& q) const {
  return {iLref_2_ * B[IR3::u],
      iLref_2_ * B[IR3::v] / (q[IR3::u] * q[IR3::u]), iLref_2_ * B[IR3::w]};
}
ddIR3 metric_cylindrical::christoffel_first_kind(const IR3& q) const {
  return {0, 0, 0, -Lref_2_ * q[IR3::u], 0, 0, 0,
      Lref_2_ * q[IR3::u], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}
ddIR3 metric_cylindrical::christoffel_second_kind(const IR3& q) const {
  return {0, 0, 0, -q[IR3::u], 0, 0, 0,
      1 / q[IR3::u], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}

}  // end namespace gyronimo
