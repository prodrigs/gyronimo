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

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/core/error.hh>
#include <gyronimo/dynamics/guiding_centre.hh>

#include <cmath>

namespace gyronimo {

//! Sets the coefficients of the equations of motion.
/*!
    Builds the guiding-centre equations of motion for a particle with
    charge-to-mass ratio `qom` (normalised to the proton's ratio) and magnetic
    moment `mu` (normalised to `Uref`/`Bref`). Providing a magnetic field is
    mandatory, the electric field is optional (it may be passed `nullptr`). If
    supplied, the optional electric field must point to the same underlying
    `metric_covariant` object as the magnetic field.
*/
guiding_centre::guiding_centre(
    double Lref, double Vref, double qom, double mu, const IR3field_c1* B,
    const IR3field* E)
    : Lref_(Lref), Vref_(Vref), qom_tilde_(qom), mu_tilde_(mu),
      magnetic_field_(B), electric_field_(E), Tref_(Lref / Vref),
      iB_time_factor_(B ? Tref_ / B->t_factor() : 0),
      iE_time_factor_(E ? Tref_ / E->t_factor() : 0),
      Oref_tilde_(
          B ? qom * codata::e / codata::m_proton * B->m_factor() * Tref_ : 1),
      iOref_tilde_(1 / Oref_tilde_) {
  if (!B) error(__func__, __FILE__, __LINE__, " no magnetic field.", 1);
  if (E && B->metric() != E->metric())
    error(__func__, __FILE__, __LINE__, " mismatched E/B coordinates.", 1);
}

//! Evaluates the time derivative `dsdt` of the dynamical state `s` at `time`.
/*!
    The dynamical system is written in general curvilinear coordinates as
    @f{equation*}{
      \frac{1}{\iota} \frac{d\tilde{q}^k}{d\tau} =
          \tilde{v}_\parallel b^k + \tilde{\Omega}^{-1} \Bigl[
              \tilde{v}_\parallel \tilde{c}^k +
                  \bigl( \mathbf{b} \times \tilde{\mathbf{d}} \bigr)^k \Bigr]
      \quad \mathrm{and} \quad
      \frac{1}{\iota} \frac{d\tilde{v}_\parallel}{d\tau} =
          -  \tilde{d}_i \Bigl( b^i + \tilde{\Omega}^{-1} \tilde{c}^i \Bigr),
    @f}
    where the vector @f$\tilde{\mathbf{d}} = \frac{1}{2} \tilde{\mu}
    \tilde{\nabla} \tilde{B} - \tilde{\mathbf{E}} + \tilde{v}_\parallel
    \partial_\tau \mathbf{b}@f$ collects all effects that cause purely
    perpendicular drifts, @f$\tilde{\mathbf{c}} = \tilde{v}_\parallel
    \tilde{\nabla} \times \mathbf{b}@f$ yields the curvature drift, and
    @f$1/\iota = 1 + \tilde{c}_\parallel/\tilde{\Omega}@f$.
*/
guiding_centre::state guiding_centre::operator()(
    const state& s, const double& time) const {
  IR3 q = this->get_position(s);
  double vpp = this->get_vpp(s);
  double jacobian = magnetic_field_->metric()->jacobian(q);
  double B_time = time * iB_time_factor_;
  IR3 covariant_b = magnetic_field_->covariant_versor(q, B_time);
  IR3 contravariant_b = magnetic_field_->contravariant_versor(q, B_time);
  auto [iO_tilde, iota, c_tilde, d_tilde] = this->dynamical_system_coefficients(
      q, vpp, B_time, jacobian, covariant_b);
  IR3 dot_X = iota *
      (vpp * contravariant_b +
       iO_tilde *
           (vpp * c_tilde +
            cross_product<contravariant>(covariant_b, d_tilde, jacobian)));
  double dot_vpp =
      -iota * inner_product(contravariant_b + iO_tilde * c_tilde, d_tilde);
  return {dot_X[IR3::u], dot_X[IR3::v], dot_X[IR3::w], dot_vpp};
}

//! Returns the sequence @f$\{1/\tilde{\Omega},\iota,\tilde{c},\tilde{d}\}@f$.
std::tuple<double, double, IR3, IR3>
guiding_centre::dynamical_system_coefficients(
    IR3& q, double vpp, double B_time, double jacobian,
    IR3& covariant_b) const {
  double inverseB = 1 / magnetic_field_->magnitude(q, B_time);
  IR3 gradB = Lref_ * magnetic_field_->del_magnitude(q, B_time);
  auto [curl_b, partial_t_b] =
      this->del_versor_b(q, B_time, inverseB, jacobian, gradB, covariant_b);
  IR3 c_tilde = vpp * curl_b;
  IR3 d_tilde = 0.5 * mu_tilde_ * gradB + vpp * partial_t_b;
  if (electric_field_) {
    double E_time = B_time * (iE_time_factor_ / iB_time_factor_);
    d_tilde -= Oref_tilde_ * electric_field_->covariant(q, E_time);
  }
  double iO_tilde = iOref_tilde_ * inverseB;
  double iota = 1 / (1 + iO_tilde * inner_product(covariant_b, c_tilde));
  return {iO_tilde, iota, c_tilde, d_tilde};
}

//! Returns @f$\{\tilde{\nabla}\times\mathbf{b},\partial_\tau\mathbf{b}\}@f$.
/*!
    Makes use of relations @f$\nabla \times \mathbf{B} = B \, \nabla \times
    \mathbf{b} + \nabla B \times \mathbf{b} @f$ and @f$\partial_t \mathbf{B} = B
    \, \partial_t \mathbf{b} + \mathbf{b} \, \partial_t B@f$ where @f$\mathbf{b}
    = \mathbf{B}/B@f$ is the versor of a vector @f$\mathbf{B}@f$.
*/
std::array<IR3, 2> guiding_centre::del_versor_b(
    IR3& q, double B_time, double inverseB, double jacobian, IR3& gradB,
    IR3& covariant_b) const {
  double partial_t_B =
      iB_time_factor_ * magnetic_field_->partial_t_magnitude(q, B_time);
  IR3 curl_b = inverseB *
      (Lref_ * magnetic_field_->curl(q, B_time) -
       cross_product<contravariant>(gradB, covariant_b, jacobian));
  IR3 partial_t_b = inverseB *
      (iB_time_factor_ * magnetic_field_->partial_t_covariant(q, B_time) -
       partial_t_B * covariant_b);
  return {curl_b, partial_t_b};
}

//! Returns the `guiding_centre::state` at a given  position, time, vpp sign.
guiding_centre::state guiding_centre::generate_state(
    const IR3& position, const double& energy_tilde,
    const guiding_centre::vpp_sign& sign, const double& time) const {
  double iLref = 1 / Lref_;
  double B_time = time * iB_time_factor_;
  double B = magnetic_field_->magnitude(position, B_time);
  double vpp = (double)(sign)*std::sqrt(energy_tilde - mu_tilde_ * B);
  return {iLref * position[0], iLref * position[1], iLref * position[2], vpp};
}

//! Returns the perpendicular energy, normalised to `Uref`.
double guiding_centre::energy_perpendicular(
    const state& s, const double& time) const {
  double B_time = time * iB_time_factor_;
  double B = magnetic_field_->magnitude(this->get_position(s), B_time);
  return mu_tilde_ * B;
}

}  // end namespace gyronimo.
