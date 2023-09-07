// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @guiding_centre.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_GUIDING_CENTRE
#define GYRONIMO_GUIDING_CENTRE

#include <gyronimo/fields/IR3field_c1.hh>

#include <array>

namespace gyronimo {

//! Guiding-centre equations of motion on a background electromagnetic field.
/*!
    Defines the equations of motion for the guiding-centre [R. Littlejohn, J.
    Plasma Phys. **29**, 111 (1983)],
    @f{gather*}{
      \biggl[ 1 + \frac{\tilde{v}_\parallel}{\tilde{\Omega}}
        \Bigl( \mathbf{b} \cdot \tilde{\nabla} \times \mathbf{b} \Bigr) \biggr]
          \frac{d\tilde{\mathbf{X}}}{d\tau} =
      \tilde{v}_\parallel \mathbf{b} + \frac{1}{\tilde{\Omega}} \biggl[
        \tilde{v}_\parallel^2 \tilde{\nabla} \times \mathbf{b} +
        \mathbf{b}\times\Bigl(
          \tilde{v}_\parallel\partial_\tau\mathbf{b} +
          \frac{1}{2}\tilde{\mu}\tilde{\nabla}\tilde{B} -
            \tilde{\mathbf{E}} \Bigr)\biggr]
      \quad \mathrm{and}\\
      \biggl[ 1 + \frac{\tilde{v}_\parallel}{\tilde{\Omega}}
        \Bigl( \mathbf{b} \cdot \tilde{\nabla} \times \mathbf{b} \Bigr) \biggr]
          \frac{d\tilde{v}_\parallel}{d\tau} = -\biggl(
        \mathbf{b} + \frac{\tilde{v}_\parallel}{\tilde{\Omega}}
          \tilde{\nabla}\times\mathbf{b} \biggr) \cdot \Bigl(
            \frac{1}{2}\tilde{\mu}\tilde{\nabla}\tilde{B} +
              \tilde{v}_\parallel\partial_\tau\mathbf{b} -
                \tilde{\mathbf{E}} \Bigr).
    @f}
    All variables are adimensional: the position of the guiding centre is
    normalised as @f$\tilde{\mathbf{X}} = \mathbf{X}/L_{ref}@f$, its parallel
    velocity as @f$\tilde{v_\parallel} = v_\parallel/V_{ref}@f$, and the time as
    @f$\tau = t/T_{ref}@f$, with the reference length `Lref` and velocity
    `Vref`=`Lref`/`Tref` being supplied to the constructor in SI units. Other
    normalisations are done internally assuming *bona-fide* electromagnetic
    fields derived from `IR3field`, with @f$\tilde{B} = B/B_{ref}@f$ whilst
    Faraday's law demands the ratio between the reference magnitudes of the
    electric and magnetic fields to match `Vref`. Other normalisations are
    @f$\tilde{\Omega} = \Omega T_{ref} = \tilde{\Omega}_{ref} \tilde{B}@f$,
    @f$\tilde{\mathbf{E}} = \tilde{\Omega}_{ref} (\mathbf{E}/E_{ref})@f$, and
    the magnetic moment is normalised to the ratio `Uref`/`Bref`, where `Uref`
    is the kinetic energy corresponding to `Vref`. Moreover, @f$\tilde{\nabla} =
    L_{ref} \nabla@f$ and @f$\mathbf{b} = \mathbf{B}/B@f$.

    The equations are implemented in a coordinate-invariant form and will work
    out-of-the-box with any coordinates defined in the `metric_covariant` object
    pointed to by the electromagnetic fields. This approach takes advantage of
    all tensor and differential calculus machinery already implemented in
    `metric_covariant`, `IR3field`, and `IR3field_c1`, which can eventually be
    specialised and optimised in further derived classes. The type
    `guiding_centre::state` implements the state of the dynamical system,
    storing the curvilinear position divided by the reference length
    (@f$\tilde{q}^\gamma = q^\gamma/L_{ref}@f$) and the normalised parallel
    velocity. Member functions are provided to convert between these types
    [i.e., `get_position(state)`, `get_vpp(state)`, `generate_state(q, v)`].
*/
class guiding_centre {
 public:
  using state = std::array<double, 4>;
  enum vpp_sign { minus = -1, plus = 1 };

  guiding_centre(
      double Lref, double Vref, double qom, double mu, const IR3field_c1* B,
      const IR3field* E);
  ~guiding_centre() {};
  state operator()(const state& s, const double& time) const;

  double Lref() const { return Lref_; };
  double Tref() const { return Tref_; };
  double Vref() const { return Vref_; };
  double mu_tilde() const { return mu_tilde_; };
  double qom_tilde() const { return qom_tilde_; };
  double Oref_tilde() const { return Oref_tilde_; };
  double get_vpp(const state& s) const { return s[3]; };
  double energy_parallel(const state& s) const;
  double energy_perpendicular(const state& s, const double& time) const;
  IR3 get_position(const state& s) const;
  state generate_state(
      const IR3& position, const double& energy_tilde, const vpp_sign& sign,
      const double& time) const;
  const IR3field* electric_field() const { return electric_field_; };
  const IR3field_c1* magnetic_field() const { return magnetic_field_; };
 private:
  const double Lref_, Vref_, qom_tilde_, mu_tilde_;
  const IR3field_c1* magnetic_field_;
  const IR3field* electric_field_;
  const double Tref_;
  const double iB_time_factor_, iE_time_factor_;
  const double Oref_tilde_, iOref_tilde_;

  std::tuple<double, double, IR3, IR3> dynamical_system_coefficients(
      IR3& q, double vpp, double B_time, double jacobian,
      IR3& covariant_b) const;
  std::array<IR3, 2> del_versor_b(
      IR3& q, double B_time, double inverseB, double jacobian, IR3& gradB,
      IR3& covariant_b) const;
};

//! Returns the parallel energy, normalised to `Uref`.
double guiding_centre::energy_parallel(const state& s) const {
  return s[3] * s[3];
}

//! Extracts curvilinear position from state.
inline IR3 guiding_centre::get_position(const state& s) const {
  return {Lref_ * s[0], Lref_ * s[1], Lref_ * s[2]};
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_GUIDING_CENTRE.
