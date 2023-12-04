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

// @morphism_cylindrical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_cylindrical.hh>

#include <cmath>

namespace gyronimo {

dIR3 morphism_cylindrical::del(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v];
  double Lref_sin = Lref_ * std::sin(phi), Lref_cos = Lref_ * std::cos(phi);
  return {Lref_cos, -r * Lref_sin, 0, Lref_sin, r * Lref_cos, 0, 0, 0, Lref_};
}
IR3 morphism_cylindrical::operator()(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v], z = q[IR3::w];
  return {Lref_ * r * std::cos(phi), Lref_ * r * std::sin(phi), Lref_ * z};
}
IR3 morphism_cylindrical::inverse(const IR3& x) const {
  double x_si = x[IR3::u], y_si = x[IR3::v], z_si = x[IR3::w];
  double r_si = std::sqrt(x_si * x_si + y_si * y_si);
  return {iLref_ * r_si, std::atan2(y_si, x_si), iLref_ * z_si};
}
ddIR3 morphism_cylindrical::ddel(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v];
  double Lref_sin = Lref_ * std::sin(phi), Lref_cos = Lref_ * std::cos(phi);
  return {
      0, -Lref_sin, 0, -r * Lref_cos, 0, 0, 0, Lref_cos,
      0, -r * Lref_sin, 0, 0, 0, 0, 0, 0, 0, 0};
}
double morphism_cylindrical::jacobian(const IR3& q) const {
  return Lref3_ * q[IR3::u];
}
dIR3 morphism_cylindrical::del_inverse(const IR3& q) const {
  double ir = 1 / q[IR3::u], phi = q[IR3::v];
  double iLref_sin = iLref_ * std::sin(phi), iLref_cos = iLref_ * std::cos(phi);
  return {
      iLref_cos, iLref_sin, 0, -iLref_sin * ir, iLref_cos * ir, 0,
      0, 0, iLref_};
}

}  // end namespace gyronimo
