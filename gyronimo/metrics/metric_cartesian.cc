// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2022 Paulo Rodrigues and Manuel Assunção.

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

// @metric_cartesian.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_cartesian.hh>

namespace gyronimo {

metric_cartesian::metric_cartesian(const morphism_cartesian* morph)
    : metric_connected(morph) {}

SM3 metric_cartesian::operator()(const IR3& r) const {
  return {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
}
SM3 metric_cartesian::inverse(const IR3& r) const {
  return {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
}
dSM3 metric_cartesian::del(const IR3& r) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
dSM3 metric_cartesian::del_inverse(const IR3& r) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
double metric_cartesian::jacobian(const IR3& r) const { return 1.0; }
IR3 metric_cartesian::del_jacobian(const IR3& r) const {
  return {0.0, 0.0, 0.0};
}
ddIR3 metric_cartesian::christoffel_first_kind(const IR3& r) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
ddIR3 metric_cartesian::christoffel_second_kind(const IR3& r) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
IR3 metric_cartesian::inertial_force(const IR3& r, const IR3& vel) const {
	return {0.0, 0.0, 0.0};
}
IR3 metric_cartesian::to_covariant(const IR3& B, const IR3& r) const {
  return B;
}
IR3 metric_cartesian::to_contravariant(const IR3& B, const IR3& r) const {
  return B;
}

}  // namespace gyronimo
