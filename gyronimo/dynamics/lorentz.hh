// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Manuel Assunção.

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

// @lorentz.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_LORENTZ
#define GYRONIMO_LORENTZ

#include <gyronimo/fields/IR3field.hh>

#include <typeinfo>

namespace gyronimo {

//! Equations of motion for a charged particle on a background electromagnetic
//! field in curvilinear coordinates.
/*!
    Defines the equations of motion for a charged particle in an electromagnetic
    field in a general set of curvilinear coordinates,
    @f{gather*}{
      \frac{dq^\alpha}{d\tau} = L_{ref} \, v^\alpha \\
      \frac{d\tilde{v}^\alpha}{d\tau} = \Omega_{ref} \left(
          \tilde{E}^\alpha + J \, \varepsilon_{ijk} \, \tilde{v}^i \,
              \tilde{B}^j \, g^{k\alpha} \right) -
                  L_{ref} \, \Gamma^\alpha_{\beta\gamma} \,
                      \tilde{v}^\beta \, \tilde{v}^\gamma
    @f}
    All variables are adimensional: the position @f$ q^\alpha @f$ of the
    particle is defined by the coordinate transformation (which can be defined
    with normalization), the time is normalized to `Tref` with @f$
    t=T_{ref}\,\tau @f$, and the velocity to `Vref`=`Lref`/`Tref` with @f$
    v^\alpha=V_{ref}\,\tilde{v}^\alpha @f$. Reference length and reference
    velocity are supplied to the constructor in SI units, other normalisations
    are done internally assuming *bona-fide* electromagnetic fields derived from
    `IR3field`, with @f$ \tilde{E} = E/E_{ref} @f$ and @f$ \tilde{B} = B/B_{ref}
    @f$ as long as Faraday's law is respected: @f$ V_{ref} = E_{ref} / B_{ref}
    @f$ where the ratio between the reference magnitudes of the electric and
    magnetic fields need to match `Vref`. Other normalisations are the reference
    frequency `Oref` with @f$ \Omega_{ref} = \frac{q_s \, B_{ref}}{m_s} \,
    T_{ref} @f$ and the reference kinetic energy `Uref` with @f$ U_{ref} =
    \frac{1}{2} m_s V_{ref}^2 @f$.

    The equations are implemented in a coordinate-invariant form and will work
    out-of-the-box with any coordinates defined in the `metric_covariant` object
    pointed to by the electromagnetic fields. This approach takes advantage of
    all tensor and differential calculus machinery already implemented in
    `metric_covariant` and `IR3field`, which can eventually be specialised and
    optimised in further derived classes. The type `lorentz::state_type`
    implements the state of the dynamical system, storing the three
    contravariant components of the particle position @f$ q^\alpha @f$ and of
    the normalised velocity @f$ \tilde{v}^\alpha @f$. Member functions are
    provided to convert between `lorentz::state_type` values, `IR3` positions,
    and `IR3` velocities [`get_position(state_type)`,
    `get_velocity(state_type)`, `generate_state(...)`].
*/
class lorentz {
 public:
  using state = std::array<double, 6>;

  lorentz(
      const double& Lref, const double& Vref, const double& qom,
      const IR3field* B, const IR3field* E);
  ~lorentz() {};

  state operator()(const state& s, const double& time) const;

  state generate_state(const IR3& q, const IR3& v) const;
  IR3 get_position(const state& s) const;
  IR3 get_velocity(const state& s) const;

  double energy_kinetic(const state& s) const;
  double energy_parallel(const state& s, const double& time) const;
  double energy_perpendicular(const state& s, const double& time) const;

  const IR3field* electric_field() const { return electric_field_; };
  const IR3field* magnetic_field() const { return magnetic_field_; };

  double Oref() const { return Oref_; };
  double Lref() const { return Lref_; };
  double Tref() const { return Tref_; };
  double Vref() const { return Vref_; };
 private:
  const double Lref_, Vref_, Tref_, Oref_;
  const IR3field *electric_field_, *magnetic_field_;
  const double iE_time_factor_, iB_time_factor_;
  const metric_covariant* metric_;

};  // end class lorentz

}  // end namespace gyronimo

#endif  // GYRONIMO_LORENTZ
