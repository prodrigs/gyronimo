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

// @metric_spherical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_spherical.hh>

#include <cmath>

namespace gyronimo {

metric_spherical::metric_spherical(const morphism_spherical* morph)
    : metric_connected(morph), Lref_(morph->Lref()), Lref2_(Lref_ * Lref_),
      Lref3_(Lref_ * Lref_ * Lref_), iLref2_(1 / Lref2_) {}
SM3 metric_spherical::operator()(const IR3& q) const {
  double factor = Lref2_ * q[IR3::u] * q[IR3::u], sinv = std::sin(q[IR3::v]);
  return {Lref2_, 0, 0, factor, 0, factor * sinv * sinv};
}
SM3 metric_spherical::inverse(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v];
  double ifactor = iLref2_ / (r * r), isin_phi = 1 / std::sin(phi);
  return {iLref2_, 0, 0, ifactor, 0, ifactor * isin_phi * isin_phi};
}
dSM3 metric_spherical::del(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v];
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double factor = 2 * Lref2_ * r;
  return {
      0, 0, 0, 0, 0, 0, 0, 0, 0, factor, 0, 0, 0, 0, 0,
      factor * sin_phi * sin_phi, factor * r * sin_phi * cos_phi, 0};
}
double metric_spherical::jacobian(const IR3& q) const {
  return Lref3_ * q[IR3::u] * q[IR3::u] * std::sin(q[IR3::v]);
}
IR3 metric_spherical::del_jacobian(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v];
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  return {Lref3_ * 2 * r * sin_phi, Lref3_ * r * r * cos_phi, 0};
}
IR3 metric_spherical::to_covariant(const IR3& B, const IR3& q) const {
  double factor = Lref2_ * q[IR3::u] * q[IR3::u];
  double sin_phi = std::sin(q[IR3::v]);
  return {
      Lref2_ * B[IR3::u], factor * B[IR3::v],
      factor * sin_phi * sin_phi * B[IR3::w]};
}
IR3 metric_spherical::to_contravariant(const IR3& B, const IR3& q) const {
  double factor = Lref2_ * q[IR3::u] * q[IR3::u], sin_phi = std::sin(q[IR3::v]);
  return {
      B[IR3::u] / Lref2_, B[IR3::v] / factor,
      B[IR3::w] / (factor * sin_phi * sin_phi)};
}
ddIR3 metric_spherical::christoffel_first_kind(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v];
  double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double factor = Lref2_ * r * sin_phi;
  double term1 =  factor * sin_phi, term2 = r * factor * cos_phi;
  return {
      0, 0, 0, -factor, 0, -term1, 0, factor, 0, 0, 0, -term2, 0, 0, term1, 0,
      term2, 0};
}
ddIR3 metric_spherical::christoffel_second_kind(const IR3& q) const {
  double r = q[IR3::u], phi = q[IR3::v];
  double ir = 1 / r, sin_phi = std::sin(phi), cos_phi = std::cos(phi);
  return {
      0, 0, 0, -r, 0, -r * sin_phi * sin_phi, 0, ir, 0, 0, 0,
      -sin_phi * cos_phi, 0, 0, ir, 0, cos_phi / sin_phi, 0};
}
IR3 metric_spherical::inertial_force(const IR3& q, const IR3& dot_q) const {
  ddIR3 gamma = christoffel_second_kind(q);
  double dot_qu = dot_q[IR3::u], dot_qv = dot_q[IR3::v], dot_qw = dot_q[IR3::w];
  return {
      -(gamma[ddIR3::uvv] * dot_qv * dot_qv +
      gamma[ddIR3::uww] * dot_qw * dot_qw),
      -(2 * gamma[ddIR3::vuv] * dot_qu * dot_qv +
      gamma[ddIR3::vww] * dot_qw * dot_qw),
      -2 * (gamma[ddIR3::wuw] * dot_qu + gamma[ddIR3::wvw] * dot_qv) * dot_qw};
}

}  // end namespace gyronimo.
