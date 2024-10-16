// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2024 Paulo Rodrigues.

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

// @aligned_frame.cc, this file is part of ::gyronimo::

#include <gyronimo/core/aligned_frame.hh>

#include <cmath>

namespace gyronimo {

aligned_frame::aligned_frame(const IR3field* field)
    : field_(field),
      metric_(dynamic_cast<const metric_connected*>(field->metric())) {
  if (!metric_) error(__func__, __FILE__, __LINE__, "no connected metric.", 1);
}

//! Frame versors {ux', uy'', uz'} at `q` after the aligning procedure.
std::array<IR3, 3> aligned_frame::versors(const IR3& q, double t) const {
  auto [b, uy_pp, ux_p] = this->gyro_versors(0, q, t);
  return {ux_p, uy_pp, b};
}

//! Gyromotion versors (Littlejohn.1983 versors {b,c,a} with b = c x a).
std::array<IR3, 3> aligned_frame::gyro_versors(
    double gyrophase, const IR3& q, double t) const {
  IR3 b = metric_->my_morphism()->from_contravariant(
      field_->contravariant_versor(q, t), q);
  double bx = b[IR3::u], by = b[IR3::v], bz = b[IR3::w];
  double r = std::hypot(bx, by);
  double cos_a = bz, sin_a = r, cos_b = by / r, sin_b = bx / r;
  double sin_phase = std::sin(gyrophase);
  double cos_phase = std::cos(gyrophase);
  IR3 vperp_versor = {
      sin_phase * cos_b - cos_phase * cos_a * sin_b,
      -sin_phase * sin_b - cos_phase * cos_a * cos_b, cos_phase * sin_a};
  IR3 rho_versor = {
      cos_phase * cos_b + sin_phase * cos_a * sin_b,
      -cos_phase * sin_b + sin_phase * cos_a * cos_b, -sin_phase * sin_a};
  return {b, vperp_versor, rho_versor};
}

//! Cartesian velocity normalised to some v_ref.
IR3 aligned_frame::velocity_from_energy_data(
    const energy_data_t& data, double gyrophase, const IR3& q, double t) const {
  double v = std::sqrt(data.energy);
  double v_parallel = v * data.pitch;
  double v_perp = v * std::sqrt(1 - std::pow(data.pitch, 2));
  auto [b, v_perp_versor, gyro_radius_versor] =
      this->gyro_versors(gyrophase, q, t);
  return v_parallel * b + data.charge_sign * v_perp * v_perp_versor;
}

//! Normalised {energy, parallel, pitch, charge_sign}.
aligned_frame::energy_data_t aligned_frame::energy_data_from_velocity(
    const IR3& v, double charge_sign, const IR3& q, double t) const {
  double energy = inner_product(v, v);  // ok because v is cartesian!
  IR3 b = metric_->my_morphism()->from_contravariant(
      field_->contravariant_versor(q, t), q);
  double pitch = inner_product(b, v) / std::sqrt(energy);
  return {energy, pitch, charge_sign};
}

}  // end namespace gyronimo.
