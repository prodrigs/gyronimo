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

// @cartesian_boris.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CARTESIAN_BORIS
#define GYRONIMO_CARTESIAN_BORIS

#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/metric_cartesian.hh>

namespace gyronimo {

//! Boris stepper, fully cartesian (position, velocity, fields).
/*!
    Advances a charged-particle's cartesian velocity by a single
    time step using the Boris algorithm [C. K. Birdsall and A. B. Langdon,
    Plasma Physics via Computer Simulation, CRC Press, 1991], which is a
    particular way of discretizing the equation of motion for the Lorentz force
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
    (B_{ref} V_{ref})@f$. Positions and velocities are cartesian and the
    supplied magnetic and electric fields are tested to check if they carry a
    `metric_cartesian` object.
*/
class cartesian_boris {
 public:
  using state = std::array<double, 6>;

  cartesian_boris(
      const double& Lref, const double& Vref, const double& qom,
      const IR3field* B, const IR3field* E);
  ~cartesian_boris() {};
  state do_step(const state& s, const double& time, const double& dt) const;
  state generate_state(const IR3& q, const IR3& v) const;
  IR3 get_position(const state& s) const;
  IR3 get_velocity(const state& s) const;

  double energy_kinetic(const state& s) const;
  double energy_parallel(const state& s, double& time) const;
  double energy_perpendicular(const state& s, double& time) const;

  state half_back_step(
      const IR3& q, const IR3& v, const double& t, const double& dt) const;

  const IR3field* electric_field() const { return electric_field_; };
  const IR3field* magnetic_field() const { return magnetic_field_; };
  const morphism* my_morphism() const { return my_morphism_; };
  double Lref() const { return Lref_; };
  double Tref() const { return Tref_; };
  double Vref() const { return Vref_; };
  double Oref() const { return Oref_; };
  double qom() const { return qom_; };
 private:
  const double Lref_, Vref_, Tref_, qom_, Oref_;
  const IR3field *electric_field_, *magnetic_field_;
  const double iE_time_factor_, iB_time_factor_;
  const metric_connected* metric_;
  const morphism* my_morphism_;
};

}  // end namespace gyronimo

#endif  // GYRONIMO_CARTESIAN_BORIS
