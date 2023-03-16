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

// @morphism_spherical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_spherical.hh>

#include <cmath>

namespace gyronimo {

IR3 morphism_spherical::operator()(const IR3& q) const {
  double r = Lref_ * q[IR3::u];
  double cn_theta = std::cos(q[IR3::v]);
  double sn_theta = std::sin(q[IR3::v]);
  double cn_phi = std::cos(q[IR3::w]);
  double sn_phi = std::sin(q[IR3::w]);
  return {r * cn_phi * sn_theta, r * sn_phi * sn_theta, r * cn_theta};
}
IR3 morphism_spherical::inverse(const IR3& x) const {
  double r_squared = x[IR3::u] * x[IR3::u] + x[IR3::v] * x[IR3::v];
  return {iLref_ * std::sqrt(r_squared + x[IR3::w] * x[IR3::w]),
      std::atan2(sqrt(r_squared), x[IR3::w]), std::atan2(x[IR3::v], x[IR3::u])};
}
dIR3 morphism_spherical::del(const IR3& q) const {
  double r = q[IR3::u];
  double cn_theta = Lref_ * std::cos(q[IR3::v]);
  double sn_theta = Lref_ * std::sin(q[IR3::v]);
  double cn_phi = std::cos(q[IR3::w]);
  double sn_phi = std::sin(q[IR3::w]);
  return {cn_phi * sn_theta, r * cn_phi * cn_theta, -r * sn_phi * sn_theta,
          sn_phi * sn_theta, r * sn_phi * cn_theta, r * cn_phi * sn_theta,
          cn_theta, -r * sn_theta, 0};
}
ddIR3 morphism_spherical::ddel(const IR3& q) const {
  double r = q[IR3::u];
  double cn_theta = Lref_ * std::cos(q[IR3::v]);
  double sn_theta = Lref_ * std::sin(q[IR3::v]);
  double cn_phi = std::cos(q[IR3::w]);
  double sn_phi = std::sin(q[IR3::w]);
  return {0, cn_phi * cn_theta, -sn_phi * sn_theta,
      -r * cn_phi * sn_theta, -r * sn_phi * cn_theta, -r * cn_phi * sn_theta,
      0, sn_phi * cn_theta, cn_phi * sn_theta,
      -r * sn_phi * sn_theta, r * cn_phi * cn_theta, -r * sn_phi * sn_theta,
      0, -sn_theta, 0,
      -r * cn_theta, 0, 0};
}
double morphism_spherical::jacobian(const IR3& q) const {
  return Lref3_ * q[IR3::u] * q[IR3::u] * std::sin(q[IR3::v]);
}
dIR3 morphism_spherical::del_inverse(const IR3& q) const {
  double ir = 1 / q[IR3::u];
  double cn_theta = iLref_ * std::cos(q[IR3::v]);
  double sn_theta = std::sin(q[IR3::v]);
  double csc_theta = iLref_ / sn_theta;
  sn_theta *= iLref_;
  double cn_phi = std::cos(q[IR3::w]);
  double sn_phi = std::sin(q[IR3::w]);
  return {cn_phi * sn_theta, sn_phi * sn_theta, cn_theta,
      ir * cn_phi * cn_theta, ir * sn_phi * cn_theta, -ir * sn_theta,
          -ir * sn_phi * csc_theta, ir * cn_phi * csc_theta, 0};
}

}  // end namespace gyronimo
