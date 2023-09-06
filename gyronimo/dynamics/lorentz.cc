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

//! Sets the coefficients of the equations of motion.
/*!
    Builds the Lorentz equations of motion for a particle with charge-to-mass
    ratio `qom` (normalised to the proton's ratio). Providing a magnetic field
    is mandatory, the electric field is optional (it may be passed `nullptr`).
    If supplied, the optional electric field must point to the same underlying
    `metric_covariant` object as the magnetic field.
*/
lorentz::lorentz(
    const double& Lref, const double& Vref, const double& qom,
    const IR3field* B, const IR3field* E)
    : Lref_(Lref), Vref_(Vref), qom_tilde_(qom), magnetic_field_(B),
      electric_field_(E), Tref_(Lref / Vref),
      iB_time_factor_(B ? Tref_ / B->t_factor() : 0),
      iE_time_factor_(E ? Tref_ / E->t_factor() : 0),
      Oref_tilde_(
          B ? qom * codata::e / codata::m_proton * B->m_factor() * Tref_ : 0),
      Eref_tilde_(
          E && B ? Oref_tilde_ * E->m_factor() / (B->m_factor() * Vref_) : 0),
      metric_(B ? B->metric() : nullptr) {
  if (!B) error(__func__, __FILE__, __LINE__, " no magnetic field.", 1);
  if (E && B->metric() != E->metric())
    error(__func__, __FILE__, __LINE__, " mismatched E/B coordinates.", 1);
}

//! Evaluates the time derivative `dsdt` of the dynamical state `s` at `time`.
/*!
    The dynamical system is written in general curvilinear coordinates as
    @f{equation*}{
      \frac{d\tilde{q}^k}{d\tau} = \tilde{v}^k
      \quad \mathrm{and} \quad
      \frac{d\tilde{v}^k}{d\tau} =
          \tilde{\Omega}_{ref}
              (\tilde{\mathbf{v}} \times \tilde{\mathbf{B}})^k +
          \tilde{E}_{ref} \tilde{E}^k
          - \tilde{\Gamma}^k_{ij} \tilde{v}^i \tilde{v}^j,
    @f}
    where @f$-\tilde{\Gamma}^k_{ij} \tilde{v}^i \tilde{v}^j@f$ are the
    inertial-force terms with @f$\tilde{\Gamma}^k_{ij} = L_{ref}
    \Gamma^k_{ij}@f$.
*/
lorentz::state lorentz::operator()(const state& s, const double& time) const {
  IR3 q = this->get_position(s), v = this->get_velocity(s);
  IR3 B = magnetic_field_->contravariant(q, iB_time_factor_ * time);
  IR3 v_cross_B = cross_product<covariant>(v, B, metric_->jacobian(q));
  IR3 dot_v = Lref_ * metric_->inertial_force(q, v) +
      Oref_tilde_ * metric_->to_contravariant(v_cross_B, q);
  if (electric_field_)
    dot_v +=
        Eref_tilde_ * electric_field_->contravariant(q, iE_time_factor_ * time);
  return {v[IR3::u],     v[IR3::v],     v[IR3::w],
          dot_v[IR3::u], dot_v[IR3::v], dot_v[IR3::w]};
}

//! Builds a state from curvilinear position `q` and normalised velocity `v`.
lorentz::state lorentz::generate_state(const IR3& q, const IR3& v) const {
  return {q[IR3::u] / Lref_, q[IR3::v] / Lref_, q[IR3::w] / Lref_,
          v[IR3::u],         v[IR3::v],         v[IR3::w]};
}

//! Returns the kinetic energy of the state, normalized to `Uref`.
double lorentz::energy_kinetic(const state& s) const {
  IR3 q = this->get_position(s), dot_tilde_q = this->get_velocity(s);
  return inner_product(dot_tilde_q, metric_->to_covariant(dot_tilde_q, q));
}

//! Returns the parallel energy of the state, normalized to `Uref`.
double lorentz::energy_parallel(const state& s, const double& time) const {
  IR3 q = this->get_position(s), dot_tilde_q = this->get_velocity(s);
  IR3 b = magnetic_field_->covariant_versor(q, iB_time_factor_ * time);
  double v_parallel = inner_product(dot_tilde_q, b);
  return v_parallel * v_parallel;
}

//! Returns the perpendicular energy of the state, normalized to `Uref`.
double lorentz::energy_perpendicular(const state& s, const double& time) const {
  IR3 q = this->get_position(s), dot_tilde_q = this->get_velocity(s);
  IR3 b = magnetic_field_->contravariant_versor(q, iB_time_factor_ * time);
  IR3 v_perpendicular =
      cross_product<covariant>(dot_tilde_q, b, metric_->jacobian(q));
  return inner_product(
      v_perpendicular, metric_->to_contravariant(v_perpendicular, q));
}

}  // end namespace gyronimo
