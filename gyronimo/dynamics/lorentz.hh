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

// @lorentz.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_LORENTZ
#define GYRONIMO_LORENTZ

#include <gyronimo/fields/IR3field.hh>

#include <typeinfo>

namespace gyronimo {

//! Lorentz-force equations of motion on a background electromagnetic field.
/*!
    Defines the equations of motion for a charged particle under the Lorentz
    force,
    @f{equation*}{
      \frac{d \tilde{\mathbf{x}}}{d\tau} = \tilde{\mathbf{v}}
      \quad \mathrm{and} \quad
      \frac{d \tilde{\mathbf{v}}}{d\tau} =
          \tilde{\Omega}_{ref} \tilde{\mathbf{v}} \times \tilde{\mathbf{B}} +
          \tilde{E}_{ref} \tilde{\mathbf{E}}.
    @f}
    All variables are adimensional: the position of the particle is normalised
    as @f$\tilde{\mathbf{x}} = \mathbf{x}/L_{ref}@f$, its velocity as
    @f$\tilde{\mathbf{v}} = \mathbf{v}/V_{ref}@f$, and the time as @f$\tau =
    t/T_{ref}@f$, with the reference length `Lref` and velocity
    `Vref`=`Lref`/`Tref` being supplied to the constructor in SI units. Other
    normalisations are done internally assuming *bona-fide* electromagnetic
    fields derived from `IR3field`, with @f$\tilde{\mathbf{B}} =
    \mathbf{B}/B_{ref}@f$ and @f$\tilde{\mathbf{E}} = \mathbf{E}/E_{ref}@f$,
    whence @f$\tilde{\Omega}_{ref} = T_{ref} \: q_s B_{ref}/m_s@f$ and
    @f$\tilde{E}_{ref} = \tilde{\Omega}_{ref} E_{ref} / (B_{ref} V_{ref})@f$.

    The equations are implemented in a coordinate-invariant form and will work
    out-of-the-box with any coordinates defined in the `metric_covariant` object
    pointed to by the electromagnetic fields. This approach takes advantage of
    all tensor and differential calculus machinery already implemented in
    `metric_covariant` and `IR3field`, which can eventually be specialised and
    optimised in further derived classes. The type `lorentz::state` implements
    the state @f$\{\tilde{q}^\gamma, \tilde{v}^\gamma \}@f$ of the dynamical
    system, storing the curvilinear position divided by the reference length
    (@f$\tilde{q}^\gamma = q^\gamma/L_{ref}@f$) and the three contravariant
    components of the normalised velocity (@f$\tilde{v}^\gamma@f$). Member
    functions are provided to convert between these types [i.e.,
    `get_position(state)`, `get_velocity(state)`, `generate_state(q, v)`].
*/
class lorentz {
 public:
  using state = std::array<double, 6>;

  lorentz(
      const double& Lref, const double& Vref, const double& qom,
      const IR3field* B, const IR3field* E);
  ~lorentz() {};
  state operator()(const state& s, const double& time) const;

  double Lref() const { return Lref_; };
  double Tref() const { return Tref_; };
  double Vref() const { return Vref_; };
  double qom_tilde() const { return qom_tilde_; };
  double Oref_tilde() const { return Oref_tilde_; };
  double Eref_tilde() const { return Eref_tilde_; };
  double energy_kinetic(const state& s) const;
  double energy_parallel(const state& s, const double& time) const;
  double energy_perpendicular(const state& s, const double& time) const;
  IR3 get_position(const state& s) const;
  IR3 get_velocity(const state& s) const;
  state generate_state(const IR3& q, const IR3& v) const;
  const IR3field* electric_field() const { return electric_field_; };
  const IR3field* magnetic_field() const { return magnetic_field_; };
 private:
  const double Lref_, Vref_, qom_tilde_;
  const IR3field *magnetic_field_, *electric_field_;
  const double Tref_;
  const double iE_time_factor_, iB_time_factor_;
  const double Oref_tilde_, Eref_tilde_;
  const metric_covariant* metric_;

};  // end class lorentz

//! Extracts the curvilinear position from a `lorentz::state`.
inline IR3 lorentz::get_position(const state& s) const {
  return {Lref_ * s[0], Lref_ * s[1], Lref_ * s[2]};
}

//! Extracts the curvilinear normalised velocity from a `lorentz::state`.
inline IR3 lorentz::get_velocity(const state& s) const {
  return {s[3], s[4], s[5]};
}

}  // end namespace gyronimo

#endif  // GYRONIMO_LORENTZ
