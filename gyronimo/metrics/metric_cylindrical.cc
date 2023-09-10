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
    : metric_connected(morph), Lref_(morph->Lref()), Lref2_(Lref_ * Lref_),
      iLref2_(1 / Lref2_), Lref3_(Lref_ * Lref_ * Lref_) {}
dSM3 metric_cylindrical::del(const IR3& q) const {
  return {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * q[IR3::u] * Lref2_, 0, 0, 0, 0, 0, 0, 0,
      0};
}
IR3 metric_cylindrical::to_covariant(const IR3& B, const IR3& q) const {
  return {
      Lref2_ * B[IR3::u], Lref2_ * B[IR3::v] * q[IR3::u] * q[IR3::u],
      Lref2_ * B[IR3::w]};
}
IR3 metric_cylindrical::to_contravariant(const IR3& B, const IR3& q) const {
  return {
      iLref2_ * B[IR3::u], iLref2_ * B[IR3::v] / (q[IR3::u] * q[IR3::u]),
      iLref2_ * B[IR3::w]};
}
ddIR3 metric_cylindrical::christoffel_first_kind(const IR3& q) const {
  return {
      0, 0, 0, -Lref2_ * q[IR3::u], 0, 0, 0, Lref2_ * q[IR3::u], 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0};
}
ddIR3 metric_cylindrical::christoffel_second_kind(const IR3& q) const {
  return {
      0, 0, 0, -q[IR3::u], 0, 0, 0, 1 / q[IR3::u], 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0};
}
IR3 metric_cylindrical::inertial_force(const IR3& q, const IR3& vel) const {
  ddIR3 gamma = christoffel_second_kind(q);
  return {
      -gamma[ddIR3::uvv] * vel[IR3::v] * vel[IR3::v],
      -2 * gamma[ddIR3::vuv] * vel[IR3::u] * vel[IR3::v], 0};
}

}  // end namespace gyronimo
