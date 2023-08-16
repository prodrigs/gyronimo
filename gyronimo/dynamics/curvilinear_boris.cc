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

// @curvilinear_boris.cc, this file is part of ::gyronimo::

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/core/error.hh>
#include <gyronimo/dynamics/boris_push.hh>
#include <gyronimo/dynamics/curvilinear_boris.hh>
#include <gyronimo/dynamics/lorentz.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

#include <boost/numeric/odeint.hpp>

#include <typeinfo>

namespace gyronimo {

//! Class constructor.
curvilinear_boris::curvilinear_boris(
    const double& Lref, const double& Vref, const double& qom,
    const IR3field* B, const IR3field* E)
    : Lref_(Lref), Vref_(Vref), Tref_(Lref / Vref), qom_(qom),
      Oref_(B ? qom * codata::e / codata::m_proton * B->m_factor() * Tref_ : 1),
      electric_field_(E), magnetic_field_(B),
      iE_time_factor_(E ? Tref_ / E->t_factor() : 1),
      iB_time_factor_(B ? Tref_ / B->t_factor() : 1),
      tildeEref_(E && B ? Oref_ * E->m_factor() / (B->m_factor() * Vref_) : 1),
      metric_(B ? dynamic_cast<const metric_connected*>(B->metric()) : nullptr),
      my_morphism_(metric_ ? metric_->my_morphism() : nullptr) {
  if (!B) error(__func__, __FILE__, __LINE__, "no magnetic field.", 1);
  if (E && B->metric() != E->metric())
    error(__func__, __FILE__, __LINE__, "mismatched E/B coordinates.", 1);
  if (!metric_)
    error(__func__, __FILE__, __LINE__, "field has no metric_connected.", 1);
}

//! Returns the update of a state by a single time step `dt`.
curvilinear_boris::state curvilinear_boris::do_step(
    const state& s, const double& time, const double& dt) const {
  double Btime = time * iB_time_factor_;

  IR3 q = this->get_position(s), v = this->get_velocity(s);
  double B = magnetic_field_->magnitude(q, Btime);
  IR3 b = my_morphism_->from_contravariant(
      magnetic_field_->contravariant_versor(q, Btime), q);
  IR3 E = electric_field_ ?
      my_morphism_->from_contravariant(
          electric_field_->contravariant(q, time * iE_time_factor_), q) :
      IR3 {0, 0, 0};

  IR3 updated_v = electric_field_ ?
      boris_push(v, Oref_, tildeEref_, E, B, b, dt) :
      boris_push(v, Oref_, B, b, dt);
  IR3 dot_q_temporary = my_morphism_->to_contravariant(updated_v, q);
  IR3 q_temporary = q + (0.5 * Lref_ * dt) * dot_q_temporary;
  IR3 v_temporary = my_morphism_->to_contravariant(updated_v, q_temporary);
  IR3 updated_q = q + (Lref_ * dt) * v_temporary;

  return this->generate_state(updated_q, updated_v);
}

//! Builds a state from curvilinear position and normalised cartesian velocity.
curvilinear_boris::state curvilinear_boris::generate_state(
    const IR3& q, const IR3& v) const {
  return {q[IR3::u], q[IR3::v], q[IR3::w], v[IR3::u], v[IR3::v], v[IR3::w]};
}

//! Extracts curvilinear position from state.
IR3 curvilinear_boris::get_position(const state& s) const {
  return {s[0], s[1], s[2]};
}

//! Extracts cartesian normalised velocity from state.
IR3 curvilinear_boris::get_velocity(const state& s) const {
  return {s[3], s[4], s[5]};
}

//! Extracts curvilinear normalised velocity from state.
IR3 curvilinear_boris::get_dot_q(const state& s) const {
  IR3 q = this->get_position(s), v = this->get_velocity(s);
  return my_morphism_->to_contravariant(v, q);
}

//! Returns the kinetic energy of the state, normalised to `Uref`.
double curvilinear_boris::energy_kinetic(const state& s) const {
  IR3 v = this->get_velocity(s);
  return inner_product(v, v);
}

//! Returns the parallel energy of the state, normalised to `Uref`.
double curvilinear_boris::energy_parallel(const state& s, double& time) const {
  IR3 q = this->get_position(s), v = this->get_velocity(s);
  IR3 b = my_morphism_->from_contravariant(
      magnetic_field_->contravariant_versor(q, time * iB_time_factor_), q);
  double v_parallel = inner_product(v, b);
  return v_parallel * v_parallel;
}

//! Returns the perpendicular energy of the state, normalised to `Uref`.
double curvilinear_boris::energy_perpendicular(
    const state& s, double& time) const {
  IR3 q = this->get_position(s), v = this->get_velocity(s);
  IR3 b = my_morphism_->from_contravariant(
      magnetic_field_->contravariant_versor(q, time * iB_time_factor_), q);
  IR3 v_perpendicular = cross_product(v, b);
  return inner_product(v_perpendicular, v_perpendicular);
}

//! Returns a state with the velocity integrated backwards by a half time step.
curvilinear_boris::state curvilinear_boris::half_back_step(
    const IR3& q, const IR3& v, const double& time, const double& dt) const {
  lorentz lo(Lref_, Vref_, qom_, magnetic_field_, electric_field_);
  auto ls = lo.generate_state(q, my_morphism_->to_contravariant(v, q));
  boost::numeric::odeint::runge_kutta4<gyronimo::lorentz::state> rk4;
  rk4.do_step(odeint_adapter(&lo), ls, time, -0.5 * dt);
  IR3 q_half_back = lo.get_position(ls), dot_q_half_back = lo.get_velocity(ls);
  IR3 v_half_back =
      my_morphism_->from_contravariant(dot_q_half_back, q_half_back);
  return this->generate_state(q, v_half_back);
}

}  // end namespace gyronimo
