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

// @classical_boris.cc, this file is part of ::gyronimo::

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/core/error.hh>
#include <gyronimo/dynamics/classical_boris.hh>
#include <gyronimo/dynamics/lorentz.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

#include <boost/numeric/odeint.hpp>

#include <cmath>

namespace gyronimo {

//! Class constructor.
/*!
    Providing a magnetic field is mandatory, the electric field is optional (it
    may be passed `nullptr`). If provided, both fields must share the same
    coordinates, i.e., the same underlying `metric_covariant` object. Such
    metric must be connected to a given morphism (and thus derive from
    `metric_connected`) in order to move between cartesian and curvilinear
    coordinates.
*/
classical_boris::classical_boris(
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
/*!
    Employs the virtual member function `morphism::translation(...)` to invert
    the specific morphism and find the updated curvilinear position.
*/
classical_boris::state classical_boris::do_step(
    const state& s, const double& time, const double& dt) const {
  IR3 q = this->get_position(s);
  IR3 updated_v = this->cartesian_velocity_update(s, time, dt);
  IR3 updated_q = my_morphism_->translation(q, (Lref_ * dt) * updated_v);
  return this->generate_state(updated_q, updated_v);
}

//! Returns the kinetic energy of the state, normalised to `Uref`.
double classical_boris::energy_kinetic(const state& s) const {
  IR3 v = this->get_velocity(s);
  return inner_product(v, v);
}

//! Returns the parallel energy of the state, normalised to `Uref`.
double classical_boris::energy_parallel(
    const state& s, const double& time) const {
  IR3 q = this->get_position(s), v = this->get_velocity(s);
  IR3 b = my_morphism_->from_contravariant(
      magnetic_field_->contravariant_versor(q, time * iB_time_factor_), q);
  double v_parallel = inner_product(v, b);
  return v_parallel * v_parallel;
}

//! Returns the perpendicular energy of the state, normalised to `Uref`.
double classical_boris::energy_perpendicular(
    const state& s, const double& time) const {
  IR3 q = this->get_position(s), v = this->get_velocity(s);
  IR3 b = my_morphism_->from_contravariant(
      magnetic_field_->contravariant_versor(q, time * iB_time_factor_), q);
  IR3 v_perpendicular = cross_product(v, b);
  return inner_product(v_perpendicular, v_perpendicular);
}

//! Returns a state with the velocity integrated backwards by a half time step.
/*!
    Integrates a Cauchy initial condition (position and velocity at the same
    instant) backwards half a time step to yield a staggered initial state.
*/
classical_boris::state classical_boris::half_back_step(
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

IR3 classical_boris::cartesian_velocity_update(
    const state& s, const double& time, const double& dt) const {
  auto [B_norm, B_versor, E] = this->get_cartesian_field_data(s, time);
  IR3 half_E_impulse = (0.5 * tildeEref_ * dt) * E;
  IR3 v_minus = this->get_velocity(s) + half_E_impulse;
  auto [T, S] = this->get_boris_rotation_coefficients(B_norm, dt);
  IR3 v_prime = v_minus + T * cross_product(v_minus, B_versor);
  IR3 v_plus = v_minus + S * cross_product(v_prime, B_versor);
  return v_plus + half_E_impulse;
}
std::array<double, 2> classical_boris::get_boris_rotation_coefficients(
    const double& B, const double& dt) const {
  double t = std::tan(0.5 * Oref_ * dt * B), s = 2 * t / (1 + t * t);
  return {t, s};
}
std::tuple<double, IR3, IR3> classical_boris::get_cartesian_field_data(
    const state& s, const double& time) const {
  IR3 q = this->get_position(s);
  IR3 E = electric_field_ ?
      my_morphism_->from_contravariant(
          electric_field_->contravariant(q, time * iE_time_factor_), q) :
      IR3 {0, 0, 0};
  IR3 B = my_morphism_->from_contravariant(
      magnetic_field_->contravariant(q, time * iB_time_factor_), q);
  double B_norm = std::sqrt(inner_product(B, B));
  return {B_norm, B / B_norm, E};
}

}  // end namespace gyronimo
