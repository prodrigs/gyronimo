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

// @cartesian_boris.cc, this file is part of ::gyronimo::

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/error.hh>
#include <gyronimo/dynamics/boris_push.hh>
#include <gyronimo/dynamics/cartesian_boris.hh>
#include <gyronimo/dynamics/lorentz.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>
#include <gyronimo/metrics/metric_cartesian.hh>

#include <boost/numeric/odeint.hpp>

namespace gyronimo {

//! Class constructor.
cartesian_boris::cartesian_boris(
    const double& Lref, const double& Vref, const double& qom,
    const IR3field* B, const IR3field* E)
    : Lref_(Lref), Vref_(Vref), Tref_(Lref / Vref), qom_(qom),
      Oref_(B ? qom * codata::e / codata::m_proton * B->m_factor() * Tref_ : 1),
      electric_field_(E), magnetic_field_(B),
      iE_time_factor_(E ? Tref_ / E->t_factor() : 1),
      iB_time_factor_(B ? Tref_ / B->t_factor() : 1),
      tildeEref_(E && B ? Oref_ * E->m_factor() / (B->m_factor() * Vref_) : 1),
      metric_(B ? dynamic_cast<const metric_cartesian*>(B->metric()) : nullptr),
      my_morphism_(metric_ ? metric_->my_morphism() : nullptr) {
  if (!B) error(__func__, __FILE__, __LINE__, "no magnetic field.", 1);
  if (E && B->metric() != E->metric())
    error(__func__, __FILE__, __LINE__, "mismatched E/B coordinates.", 1);
  if (!metric_)
    error(__func__, __FILE__, __LINE__, "field has no metric_cartesian.", 1);
}

//! Performs a time step `dt` update to the state `in` and returns the result.
cartesian_boris::state cartesian_boris::do_step(
    const state& s, const double& time, const double& dt) const {
  double Btime = time * iB_time_factor_;

  IR3 x = this->get_position(s), v = this->get_velocity(s);
  double B = magnetic_field_->magnitude(x, Btime);
  IR3 b = magnetic_field_->contravariant_versor(x, Btime);
  IR3 E = electric_field_ ?
      electric_field_->contravariant(x, time * iE_time_factor_) :
      IR3 {0, 0, 0};

  IR3 updated_v = electric_field_ ?
      boris_push(v, Oref_, tildeEref_, E, B, b, dt) :
      boris_push(v, Oref_, B, b, dt);
  IR3 updated_x = x + (Lref_ * dt) * updated_v;

  return this->generate_state(updated_x, updated_v);
}

//! Returns a `cartesian_boris::state` from SI position and normalised velocity.
cartesian_boris::state cartesian_boris::generate_state(
    const IR3& x, const IR3& v) const {
  return {x[IR3::u], x[IR3::v], x[IR3::w], v[IR3::u], v[IR3::v], v[IR3::w]};
}

//! Returns the vector position of the state in SI units.
IR3 cartesian_boris::get_position(const state& s) const {
  return {s[0], s[1], s[2]};
}

//! Returns the vector velocity of the state (normalised to `Lref/Tref`).
IR3 cartesian_boris::get_velocity(const state& s) const {
  return {s[3], s[4], s[5]};
}

//! Returns the kinetic energy of the state, normalized to `Uref`.
double cartesian_boris::energy_kinetic(const state& s) const {
  IR3 v = this->get_velocity(s);
  return inner_product(v, v);
}

//! Returns the parallel energy of the state, normalized to `Uref`.
double cartesian_boris::energy_parallel(const state& s, double& time) const {
  IR3 x = this->get_position(s), v = this->get_velocity(s);
  IR3 b = magnetic_field_->contravariant_versor(x, time * iB_time_factor_);
  double v_parallel = inner_product(v, b);
  return v_parallel * v_parallel;
}

//! Returns the perpendicular energy of the state, normalized to `Uref`.
double cartesian_boris::energy_perpendicular(
    const state& s, double& time) const {
  IR3 x = this->get_position(s), v = this->get_velocity(s);
  IR3 b = magnetic_field_->contravariant_versor(x, time * iB_time_factor_);
  IR3 v_perpendicular = cross_product(v, b);
  return inner_product(v_perpendicular, v_perpendicular);
}

//! Returns a state with the velocity integrated backwards by a half time step.
cartesian_boris::state cartesian_boris::half_back_step(
    const IR3& x, const IR3& v, const double& time, const double& dt) const {
  lorentz lo(Lref_, Vref_, qom_, magnetic_field_, electric_field_);
  auto ls = lo.generate_state(x, v);
  boost::numeric::odeint::runge_kutta4<gyronimo::lorentz::state> rk4;
  rk4.do_step(odeint_adapter(&lo), ls, time, -0.5 * dt);
  IR3 v_half_back = lo.get_velocity(ls);
  return this->generate_state(x, v_half_back);
}

}  // end namespace gyronimo
