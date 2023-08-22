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

// @curvilinear_boris.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CURVILINEAR_BORIS
#define GYRONIMO_CURVILINEAR_BORIS

#include <gyronimo/dynamics/classical_boris.hh>

namespace gyronimo {

//! Classical Boris-like stepper with alternative curvilinear-position advance.
/*!
    Advances a charged-particle's cartesian velocity by a single time step using
    the same algorithm implemented in `classical_boris` to discretise the
    equation of motion for the Lorentz force
    @f{equation*}{
      \frac{d\tilde{\mathbf{v}}}{d\tau}
          = \tilde{\Omega}_{ref} \tilde{\mathbf{v}} \times \tilde{\mathbf{B}} +
            \tilde{E}_{ref} \tilde{\mathbf{E}}.
    @f}
    However, the curvilinear position is here advanced using the technique
    described in G. Delzanno *et al*., IEEE Trans. Plasma Sci. **41**, 3577
    (2013): the Runge/Kutta-2 midpoint method provides a temporary estimate
    @f{equation*}{
      q^k_{\tau + \Delta\tau/2} = q^k_\tau + \frac{1}{2} \Delta\tau \:
          \underbrace{
          \tilde{\mathbf{v}}_{\tau + \Delta\tau/2} \cdot \mathbf{e}^k_\tau
          }_{\dot{q}^k_\star},
    @f}
    which is used to evaluate the tangent-space basis @f$\mathbf{e}^k(q^k)@f$ at
    the velocity instant in order to produce
    @f{equation*}{
      q^k_{\tau + \Delta\tau} = q^k_\tau + \Delta\tau \:
          \underbrace{
          \tilde{\mathbf{v}}_{\tau + \Delta\tau/2} \cdot
          \mathbf{e}^k_{\tau + \Delta\tau/2}}_{\dot{q}^k_{\tau + \Delta\tau/2}}.
    @f}
    This approach avoids the numerical cost usually required to invert the
    morphism associated with the specific coordinates being used. Notice that
    such cost may eventually be neglegible if the specific morphism is inverted
    analytically.
*/
class curvilinear_boris {
 public:
  using state = classical_boris::state;

  curvilinear_boris(
      const double& Lref, const double& Vref, const double& qom,
      const IR3field* B, const IR3field* E)
      : classical_boris_(Lref, Vref, qom, B, E) {};
  ~curvilinear_boris() {};
  state do_step(const state& s, const double& time, const double& dt) const;

  double Lref() const { return classical_boris_.Lref(); };
  double Tref() const { return classical_boris_.Tref(); };
  double Vref() const { return classical_boris_.Vref(); };
  double Oref() const { return classical_boris_.Oref(); };
  double qom() const { return classical_boris_.qom(); };
  double energy_kinetic(const state& s) const;
  double energy_parallel(const state& s, const double& time) const;
  double energy_perpendicular(const state& s, const double& time) const;
  IR3 get_position(const state& s) const;
  IR3 get_velocity(const state& s) const;
  IR3 get_dot_q(const state& s) const;
  state generate_state(const IR3& q, const IR3& v) const;
  const IR3field* electric_field() const;
  const IR3field* magnetic_field() const;
  const morphism* my_morphism() const;

  IR3 cartesian_velocity_update(
      const state& s, const double& t, const double& dt) const;
  state half_back_step(
      const IR3& q, const IR3& v, const double& t, const double& dt) const;
 private:
  const classical_boris classical_boris_;
};

inline double curvilinear_boris::energy_kinetic(const state& s) const {
  return classical_boris_.energy_kinetic(s);
}
inline double curvilinear_boris::energy_parallel(
    const state& s, const double& time) const {
  return classical_boris_.energy_parallel(s, time);
}
inline double curvilinear_boris::energy_perpendicular(
    const state& s, const double& time) const {
  return classical_boris_.energy_perpendicular(s, time);
}
inline IR3 curvilinear_boris::get_position(const state& s) const {
  return classical_boris_.get_position(s);
}
inline IR3 curvilinear_boris::get_velocity(const state& s) const {
  return classical_boris_.get_velocity(s);
}
inline IR3 curvilinear_boris::get_dot_q(const state& s) const {
  return classical_boris_.get_dot_q(s);
};
inline curvilinear_boris::state curvilinear_boris::generate_state(
    const IR3& q, const IR3& v) const {
  return classical_boris_.generate_state(q, v);
}
inline const IR3field* curvilinear_boris::electric_field() const {
  return classical_boris_.electric_field();
}
inline const IR3field* curvilinear_boris::magnetic_field() const {
  return classical_boris_.magnetic_field();
}
inline const morphism* curvilinear_boris::my_morphism() const {
  return classical_boris_.my_morphism();
}
inline IR3 curvilinear_boris::cartesian_velocity_update(
    const state& s, const double& t, const double& dt) const {
  return classical_boris_.cartesian_velocity_update(s, t, dt);
}
inline curvilinear_boris::state curvilinear_boris::half_back_step(
    const IR3& q, const IR3& v, const double& t, const double& dt) const {
  return classical_boris_.half_back_step(q, v, t, dt);
}

}  // end namespace gyronimo

#endif  // GYRONIMO_CURVILINEAR_BORIS
