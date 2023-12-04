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

// @morphism_polar_torus.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_polar_torus.hh>

#include <cmath>

namespace gyronimo {

morphism_polar_torus::morphism_polar_torus(
    const double minor_radius, double major_radius)
    : minor_radius_(minor_radius), major_radius_(major_radius),
      iaspect_ratio_(minor_radius_ / major_radius_),
      volume_factor_(minor_radius * minor_radius * major_radius),
      iminor_radius_(1.0 / minor_radius) {}
IR3 morphism_polar_torus::operator()(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v], phi = q[IR3::w];
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double R = major_radius_ * (1.0 + iaspect_ratio_ * r * cos_theta);
  return {R * cos_phi, -R * sin_phi, minor_radius_ * r * sin_theta};
}
IR3 morphism_polar_torus::inverse(const IR3& x) const {
  double x_si = x[IR3::u], y_si = x[IR3::v], z_si = x[IR3::w];
  double R = std::sqrt(x_si * x_si + y_si * y_si);
  double deltaR = R - major_radius_;
  return {
      iminor_radius_ * std::sqrt(z_si * z_si + deltaR * deltaR),
      std::atan2(z_si, deltaR), std::atan2(-y_si, x_si)};
}
dIR3 morphism_polar_torus::del(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v], phi = q[IR3::w];
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  double R = major_radius_ * (1.0 + iaspect_ratio_ * r * cos_theta);
  double a_cos_theta = minor_radius_ * cos_theta;
  double ar_cos_theta = r * a_cos_theta;
  double a_sin_theta = minor_radius_ * sin_theta;
  double ar_sin_theta = r * a_sin_theta;
  return {
      a_cos_theta * cos_phi, -ar_sin_theta * cos_phi, -R * sin_phi,
      -a_cos_theta * sin_phi, ar_sin_theta * sin_phi, -R * cos_phi,
      a_sin_theta, ar_cos_theta, 0};
}
ddIR3 morphism_polar_torus::ddel(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v], phi = q[IR3::w];
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  double R = major_radius_ * (1.0 + iaspect_ratio_ * r * cos_theta);
  double a_cos_theta = minor_radius_ * cos_theta;
  double ar_cos_theta = r * a_cos_theta;
  double a_sin_theta = minor_radius_ * sin_theta;
  double ar_sin_theta = r * a_sin_theta;
  return {
      0, -a_sin_theta * cos_phi, -a_cos_theta * sin_phi,
      -ar_cos_theta * cos_phi, ar_sin_theta * sin_phi, -R * cos_phi, 0,
      a_sin_theta * sin_phi, -a_cos_theta * cos_phi, ar_cos_theta * sin_phi,
      ar_sin_theta * cos_phi, R * sin_phi, 0, a_cos_theta, 0, -ar_sin_theta,
      0, 0};
}
double morphism_polar_torus::jacobian(const IR3& q) const {
  double r = q[IR3::u], cos_theta = std::cos(q[IR3::v]);
  double Rfactor = (1.0 + iaspect_ratio_ * r * cos_theta);
  return volume_factor_ * r * Rfactor;
}
dIR3 morphism_polar_torus::del_inverse(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v], phi = q[IR3::w], ir = 1 / r;
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  double iR = 1.0 / (major_radius_ * (1.0 + iaspect_ratio_ * r * cos_theta));
  double ia_cos_theta = iminor_radius_ * cos_theta;
  double iar_cos_theta = ir * ia_cos_theta;
  double ia_sin_theta = iminor_radius_ * sin_theta;
  double iar_sin_theta = ir * ia_sin_theta;
  return {
      ia_cos_theta * cos_phi, -ia_cos_theta * sin_phi, ia_sin_theta,
      -iar_sin_theta * cos_phi, iar_sin_theta * sin_phi, iar_cos_theta,
      -iR * sin_phi, -iR * cos_phi, 0};
}

}  // end namespace gyronimo
