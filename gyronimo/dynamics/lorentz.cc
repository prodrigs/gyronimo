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

// @lorentz.cc, this file is part of ::gyronimo::

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/core/error.hh>
#include <gyronimo/dynamics/lorentz.hh>

#include <cmath>

namespace gyronimo {

//! Class constructor.
/*!
    Builds the Lorentz equations of motion for a particle with charge-to-mass
    ratio `qom` (normalised to the proton's ratio). If supplied (defaults to no
    field), the optional electric field must point to the same underlying
    `metric_covariant` object as the magnetic field.
*/
lorentz::lorentz(
    const double& Lref, const double& Vref, const double& qom,
    const IR3field* B, const IR3field* E)
    : Lref_(Lref), Vref_(Vref), Tref_(Lref / Vref),
      Oref_(
          B ? qom * codata::e / codata::m_proton * B->m_factor() * Tref_ : 1.0),
      electric_field_(E), magnetic_field_(B),
      iE_time_factor_(E ? Tref_ / E->t_factor() : 1.0),
      iB_time_factor_(B ? Tref_ / B->t_factor() : 1.0),
      metric_(B ? B->metric() : nullptr) {
  if (!B) error(__func__, __FILE__, __LINE__, " no magnetic field.", 1);
  if (E && B->metric() != E->metric())
    error(__func__, __FILE__, __LINE__, " mismatched E/B coordinates.", 1);
}

//! Evaluates the time derivative `dsdt` of the dynamical state `s` at time `t`.
lorentz::state lorentz::operator()(const state& s, const double& time) const {
  IR3 q = {s[0], s[1], s[2]};
  IR3 v = {s[3], s[4], s[5]};

  IR3 dot_q = Lref_ * v;
  IR3 dot_v = Lref_ * metric_->inertial_force(q, v);

  IR3 B = magnetic_field_->contravariant(q, iB_time_factor_ * time);
  IR3 vxB = cross_product<covariant>(v, B, metric_->jacobian(q));
  dot_v += Oref_ * metric_->to_contravariant(vxB, q);
  if (electric_field_)
    dot_v += Oref_ * electric_field_->contravariant(q, iE_time_factor_ * time);

  return {dot_q[IR3::u], dot_q[IR3::v], dot_q[IR3::w],
          dot_v[IR3::u], dot_v[IR3::v], dot_v[IR3::w]};
}

//! Generates a `lorentz::state` from curvilinear position `q` and velocity `v`.
lorentz::state lorentz::generate_state(const IR3& q, const IR3& v) const {
  return {q[IR3::u], q[IR3::v], q[IR3::w], v[IR3::u], v[IR3::v], v[IR3::w]};
}

//! Extracts the position from a `lorentz::state`.
IR3 lorentz::get_position(const state& s) const { return {s[0], s[1], s[2]}; }

//! Extracts the velocity from a `lorentz::state`, normalized to `Vref`.
IR3 lorentz::get_velocity(const state& s) const { return {s[3], s[4], s[5]}; }

//! Returns the kinetic energy of the state, normalized to `Uref`.
double lorentz::energy_kinetic(const state& s) const {
  IR3 q = {s[0], s[1], s[2]};
  IR3 dot_q = {s[3], s[4], s[5]};
  return inner_product(dot_q, metric_->to_covariant(dot_q, q));
}

//! Returns the parallel energy of the state, normalized to `Uref`.
double lorentz::energy_parallel(const state& s, const double& time) const {
  IR3 q = {s[0], s[1], s[2]};
  IR3 dot_q = {s[3], s[4], s[5]};
  IR3 b = magnetic_field_->covariant_versor(q, iB_time_factor_ * time);
  double vpar = inner_product(dot_q, b);
  return vpar * vpar;
}

//! Returns the perpendicular energy of the state, normalized to `Uref`.
double lorentz::energy_perpendicular(const state& s, const double& time) const {
  IR3 q = {s[0], s[1], s[2]};
  IR3 dot_q = {s[3], s[4], s[5]};
  IR3 b = magnetic_field_->contravariant_versor(q, iB_time_factor_ * time);
  IR3 vperp = cross_product<covariant>(dot_q, b, metric_->jacobian(q));
  return inner_product(vperp, metric_->to_contravariant(vperp, q));
}

}  // end namespace gyronimo
