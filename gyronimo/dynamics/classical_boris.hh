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

// @classical_boris.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CLASSICAL_BORIS
#define GYRONIMO_CLASSICAL_BORIS

#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Classical Boris-like stepper, cartesian velocity and curvilinear position.
/*!
    Advances a charged-particle's cartesian velocity by a single time step using
    the Boris algorithm [C. K. Birdsall and A. B. Langdon, Plasma Physics via
    Computer Simulation, CRC Press, 1991], which is a particular way of
    discretizing the equation of motion for the Lorentz force
    @f{equation*}{
      \frac{d\tilde{\mathbf{v}}}{d\tau}
          = \tilde{\Omega}_{ref} \tilde{\mathbf{v}} \times \tilde{\mathbf{B}} +
            \tilde{E}_{ref} \tilde{\mathbf{E}}
    @f}
    Here, @f$\tau = t/T_{ref}@f$, @f$\tilde{\mathbf{v}} = \mathbf{v}/V_{ref}@f$
    @f$\tilde{\mathbf{B}} = \mathbf{B}/B_{ref}@f$, and @f$\tilde{\mathbf{E}} =
    \mathbf{E}/E_{ref}@f$ are the time, cartesian velocity, magnetic, and
    electric fields normalised, respectively, to `Tref`, `Vref`, `Bref`, and
    `Eref` (all in SI units), whereas @f$\tilde{\Omega}_{ref} = T_{ref} \: q_s
    B_{ref}/m_s@f$ and @f$\tilde{E}_{ref} = \tilde{\Omega}_{ref} E_{ref} /
    (B_{ref} V_{ref})@f$. The curvilinear position is advanced from
    @f$q^k_\tau@f$ to @f$q^k_{\tau + \Delta\tau}@f$ by solving the (possibly
    nonlinear) system
    @f{equation*}{
      V_{ref} T_{ref} \: \Delta\tau \: \tilde{\mathbf{v}}_{\tau + \Delta\tau/2}
          = \mathcal{M}\big(q^k_{\tau + \Delta\tau}\big) -
              \mathcal{M}\big(q^k_\tau\big)
    @f}
    where @f$\mathcal{M}(q^k)@f$ is the morphism from the coordinates @f$q^k@f$
    into the cartesian space carried by the specific magnetic and electric field
    objects supplied to the constructor. Notice that the conventional Boris
    stepper in cartesian coordinates is recovered if @$f\mathcal{M}@$f is the
    identity (e.g., `morphism_cartesian`).
*/
class classical_boris {
 public:
  using state = std::array<double, 6>;

  classical_boris(
      const double& Lref, const double& Vref, const double& qom,
      const IR3field* B, const IR3field* E);
  ~classical_boris() {};
  state do_step(const state& s, const double& time, const double& dt) const;

  double Lref() const { return Lref_; };
  double Tref() const { return Tref_; };
  double Vref() const { return Vref_; };
  double Oref() const { return Oref_; };
  double qom() const { return qom_; };
  double energy_kinetic(const state& s) const;
  double energy_parallel(const state& s, const double& time) const;
  double energy_perpendicular(const state& s, const double& time) const;
  IR3 get_position(const state& s) const { return {s[0], s[1], s[2]}; };
  IR3 get_velocity(const state& s) const { return {s[3], s[4], s[5]}; };
  IR3 get_dot_q(const state& s) const;
  state generate_state(const IR3& q, const IR3& v) const;
  const IR3field* electric_field() const { return electric_field_; };
  const IR3field* magnetic_field() const { return magnetic_field_; };
  const morphism* my_morphism() const { return my_morphism_; };

  IR3 cartesian_velocity_update(
      const state& s, const double& time, const double& dt) const;
  state half_back_step(
      const IR3& q, const IR3& v, const double& t, const double& dt) const;
 private:
  const double Lref_, Vref_, Tref_, qom_, Oref_;
  const IR3field *electric_field_, *magnetic_field_;
  const double iE_time_factor_, iB_time_factor_, tildeEref_;
  const metric_connected* metric_;
  const morphism* my_morphism_;

  std::array<double, 2> get_boris_rotation_coefficients(
      const double& B, const double& dt) const;
  std::tuple<double, IR3, IR3> get_cartesian_field_data(
      const state& s, const double& time) const;
};

//! Extracts curvilinear normalised velocity from state.
inline IR3 classical_boris::get_dot_q(const state& s) const {
  IR3 q = this->get_position(s), v = this->get_velocity(s);
  return my_morphism_->to_contravariant(v, q);
}

//! Builds a state from curvilinear position and normalised cartesian velocity.
inline classical_boris::state classical_boris::generate_state(
    const IR3& q, const IR3& v) const {
  return {q[IR3::u], q[IR3::v], q[IR3::w], v[IR3::u], v[IR3::v], v[IR3::w]};
}

}  // end namespace gyronimo

#endif  // GYRONIMO_CLASSICAL_BORIS
