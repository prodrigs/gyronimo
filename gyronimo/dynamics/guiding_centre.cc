// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues.

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

// @guiding_centre.cc, this file is part of ::gyronimo::

#include <cmath>
#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/dynamics/guiding_centre.hh>

namespace gyronimo {

//! Class constructor.
/*!
    Builds the guiding-centre equations for a particle with charge-to-mass ratio
    `qom` (normalised to the proton's ratio) and magnetic moment `mu`
    (normalised to `Uref`/`Bref`). The magnetic field must be continuously
    differentiable whilst the electric field needs only to be continuous (it is
    also optional, defaults to no field).
*/
guiding_centre::guiding_centre(
    double Lref, double Vref,
    double qom, double mu, const IR3field_c1* B, const IR3field* E)
    : qom_tilde_(qom), mu_tilde_(mu),
      electric_field_(E), magnetic_field_(B),
      Lref_(Lref), Vref_(Vref), Tref_(Lref/Vref),
      Bfield_time_factor_(Lref/(Vref*B->t_factor())),
      Efield_time_factor_(E ? Lref/(Vref*E->t_factor()) : 0) {
  iOref_tilde_ = 1.0/(qom*codata::e/codata::m_proton*B->m_factor()*Tref_);
  Eref_tilde_ = (E ? E->m_factor() : 0)/(Vref_*B->m_factor()*iOref_tilde_);
}

//! Evaluates the time derivative `dxdt` of the dynamical state `x` at time `t`.
guiding_centre::state guiding_centre::operator()(
    const state& x, double t) const {
  double Btime = t*Bfield_time_factor_;
  double Etime = t*Efield_time_factor_;
  double vpp = this->get_vpp(x);
  IR3 position = this->get_position(x);

  double inverseB = 1.0/magnetic_field_->magnitude(position, Btime);
  double jacobian = magnetic_field_->metric()->jacobian(position);
  double partial_t_B = Bfield_time_factor_*
      magnetic_field_->partial_t_magnitude(position, Btime);
  IR3 gradB = Lref_*magnetic_field_->del_magnitude(position, Btime);
  IR3 covariant_b = magnetic_field_->covariant_versor(position, Btime);
  IR3 contravariant_b = magnetic_field_->contravariant_versor(position, Btime);

// Makes use of relations curl(\vec{B}) = B curl(\vec{b}) + grad(B)x\vec{b}
// and \partial_t(\vec{B}) = B \partial_t(\vec{b}) + \vec{b}\partial_t B.
  IR3 curlb = inverseB*(
      Lref_*magnetic_field_->curl(position, Btime) -
      cross_product(gradB, covariant_b, jacobian));
  IR3 partial_t_b = inverseB*(
      Bfield_time_factor_*magnetic_field_->partial_t_covariant(
          position, Btime) -
              partial_t_B*magnetic_field_->covariant_versor(position, Btime));

  double inverse_omega_tilde = iOref_tilde_*inverseB;
  double prefactor = 1.0/(
      1.0 + inverse_omega_tilde*vpp*inner_product(covariant_b, curlb));
  IR3 drifter = 0.5*mu_tilde_*gradB + vpp*partial_t_b;
  if (electric_field_)
      drifter -= Eref_tilde_*electric_field_->covariant(position, Etime);

  IR3 dot_X = prefactor*(vpp*contravariant_b +
      inverse_omega_tilde*(
          vpp*vpp*curlb + cross_product(covariant_b, drifter, jacobian)));
  double dot_vpp = -prefactor*inner_product(
      contravariant_b + inverse_omega_tilde*vpp*curlb, drifter);
  return {dot_X[IR3::u], dot_X[IR3::v], dot_X[IR3::w], dot_vpp};
}

double guiding_centre::get_vpp(const state& s) const {
  return s[3];
}

IR3 guiding_centre::get_position(const state& s) const {
  return {Lref_*s[0], Lref_*s[1], Lref_*s[2]};
}

//! Returns the `guiding_centre::state` at a given  position, time, vpp sign.
guiding_centre::state guiding_centre::generate_state(
    const IR3& position, double energy_tilde,
    guiding_centre::vpp_sign sign, double time) const {
  double iLref = 1.0/Lref_;
  double Btime = time*Bfield_time_factor_;
  double B = magnetic_field_->magnitude(position, Btime);
  double vpp = (double)(sign)*std::sqrt(energy_tilde - mu_tilde_*B);
  return {iLref*position[0], iLref*position[1], iLref*position[2], vpp};
}

//! Returns the parallel energy, normalised to `Uref`.
double guiding_centre::energy_parallel(const state &s) const {
  return s[3]*s[3];
}

//! Returns the perpendicular energy, normalised to `Uref`.
double guiding_centre::energy_perpendicular(const state &s, double time) const {
  double Btime = time*Bfield_time_factor_;
  double B = magnetic_field_->magnitude(this->get_position(s), Btime);
  return mu_tilde_*B;
}

} // end namespace gyronimo.
