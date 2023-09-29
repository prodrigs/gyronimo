// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 20222023 Manuel Assunção and Paulo Rodrigues.

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

// @morphism_spherical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_spherical.hh>

#include <cmath>

namespace gyronimo {

IR3 morphism_spherical::operator()(const IR3& q) const {
  double r_si = Lref_ * q[IR3::u], phi = q[IR3::v], theta = q[IR3::w];
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  return {
      r_si * cos_theta * sin_phi, r_si * sin_theta * sin_phi, r_si * cos_phi};
}
IR3 morphism_spherical::inverse(const IR3& x) const {
  double x_si = x[IR3::u], y_si = x[IR3::v], z_si = x[IR3::w];
  double r_sin_phi_squared = x_si * x_si + y_si * y_si;
  return {
      iLref_ * std::sqrt(r_sin_phi_squared + z_si * z_si),
      std::atan2(sqrt(r_sin_phi_squared), z_si), std::atan2(y_si, x_si)};
}
dIR3 morphism_spherical::del(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v], theta = q[IR3::w];
  double cos_phi = Lref_ * std::cos(phi), sin_phi = Lref_ * std::sin(phi);
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  return {
      cos_theta * sin_phi, r * cos_theta * cos_phi, -r * sin_theta * sin_phi,
      sin_theta * sin_phi, r * sin_theta * cos_phi, r * cos_theta * sin_phi,
      cos_phi, -r * sin_phi, 0};
}
ddIR3 morphism_spherical::ddel(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v], theta = q[IR3::w];
  double cos_phi = Lref_ * std::cos(phi), sin_phi = Lref_ * std::sin(phi);
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  return {
      0, cos_theta * cos_phi, -sin_theta * sin_phi, -r * cos_theta * sin_phi,
      -r * sin_theta * cos_phi, -r * cos_theta * sin_phi, 0,
      sin_theta * cos_phi, cos_theta * sin_phi, -r * sin_theta * sin_phi,
      r * cos_theta * cos_phi, -r * sin_theta * sin_phi, 0, -sin_phi, 0,
      -r * cos_phi, 0, 0};
}
dIR3 morphism_spherical::del_inverse(const IR3& q) const {
  double ir = 1 / q[IR3::u], phi = q[IR3::v], theta = q[IR3::w];
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  double cos_phi = iLref_ * std::cos(phi), sin_phi = std::sin(phi),
         csc_phi = iLref_ / sin_phi;
  sin_phi *= iLref_;
  return {
      cos_theta * sin_phi, sin_theta * sin_phi, cos_phi,
      ir * cos_theta * cos_phi, ir * sin_theta * cos_phi, -ir * sin_phi,
      -ir * sin_theta * csc_phi, ir * cos_theta * csc_phi, 0};
}
double morphism_spherical::jacobian(const IR3& q) const {
  return Lref3_ * q[IR3::u] * q[IR3::u] * std::sin(q[IR3::v]);
}

}  // end namespace gyronimo
