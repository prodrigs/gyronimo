// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues and Manuel Assunção.

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

// @metric_polar_torus.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_polar_torus.hh>

#include <cmath>

namespace gyronimo {

metric_polar_torus::metric_polar_torus(const morphism_polar_torus* morph)
    : metric_connected(morph), minor_radius_(morph->minor_radius()),
      major_radius_(morph->major_radius()),
      minor_radius_squared_(minor_radius_ * minor_radius_),
      iminor_radius_squared_(1 / minor_radius_squared_),
      major_radius_squared_(major_radius_ * major_radius_),
      iaspect_ratio_(morph->iaspect_ratio()) {}
SM3 metric_polar_torus::operator()(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v];
  double R = major_radius_ * (1 + iaspect_ratio_ * r * std::cos(theta));
  return {minor_radius_squared_, 0, 0, minor_radius_squared_ * r * r, 0, R * R};
}
SM3 metric_polar_torus::inverse(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v];
  double R = major_radius_ * (1 + iaspect_ratio_ * r * std::cos(theta));
  return {
      iminor_radius_squared_, 0, 0, iminor_radius_squared_ / (r * r), 0,
      1 / (R * R)};
}
dSM3 metric_polar_torus::del(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v];
  double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  double R = major_radius_ * (1 + iaspect_ratio_ * r * cos_theta);
  double factor = 2 * R * minor_radius_;
  return {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 2 * r * minor_radius_squared_, 0, 0, 0, 0, 0,
      factor * cos_theta, -factor * r * sin_theta, 0};
}
double metric_polar_torus::jacobian(const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v];
  return minor_radius_squared_ * major_radius_ * r *
      (1 + iaspect_ratio_ * r * std::cos(theta));
}
IR3 metric_polar_torus::to_covariant(const IR3& B, const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v];
  double R = major_radius_ * (1 + iaspect_ratio_ * r * std::cos(theta));
  return {
      minor_radius_squared_ * B[IR3::u],
      minor_radius_squared_ * r * r * B[IR3::v], R * R * B[IR3::w]};
}
IR3 metric_polar_torus::to_contravariant(const IR3& B, const IR3& q) const {
  double r = q[IR3::u], theta = q[IR3::v];
  double R = major_radius_ * (1 + iaspect_ratio_ * r * std::cos(theta));
  return {
      B[IR3::u] / minor_radius_squared_,
      B[IR3::v] / (minor_radius_squared_ * r * r), B[IR3::w] / (R * R)};
}

}  // end namespace gyronimo.
