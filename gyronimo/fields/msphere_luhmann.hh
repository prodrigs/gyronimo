// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues.

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

// @msphere_luhmann.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MSPHERE_LUHMANN
#define GYRONIMO_MSPHERE_LUHMANN

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/metrics/metric_spherical.hh>

namespace gyronimo {

//! Simple analytical model of the Earth magnetosphere.
/*!
    Implements a linear combination of a dipole and a current-sheet field [J.
    Luhmann and L. Friesen, J. Geophys. Res. **84**, 4405 (1979)] using
    spherical coordinates where the geocentric distance is normalised to
    `earth_radius` (i.e., @f$\tilde{r} = r/R_E@f$) and the azimuthal angle is
    measured from noon. By default, the field is normalised to
    `earth_surface_avg_field`. The magnitudes of the dipole and current-sheet
    fields may be provided in units of @f$\mathrm{Gauss} \times R_E^3@f$ and
    mGauss respectively, the default (and recommended) values being
    `dipole_factor = 0.31` and `csheet_factor = 0.15`. The current sheet at the
    equatorial plane is smoothed by a tanh function with a `smooth_factor`
    (i.e., @f$\delta@f$) provided in @f$R_E@f$ units. Overall, one has
    @f{equation*}{
      \mathbf{B} =
        d \nabla \biggl(\frac{\cos \theta}{\tilde{r}^2} \biggr) +
        c \tanh \biggl(\frac{\tilde{r} \cos \theta}{\delta} \biggr)
          \mathbf{u}_x
    @f}
    where the versor @f$\mathbf{u}_x@f$ is directed towards the Sun.
*/
class msphere_luhmann : public IR3field_c1 {
 public:
  static constexpr double earth_radius = 6378137.0;  // [SI]
  static constexpr double earth_surface_avg_field = 0.5e-04;  // [SI]
  msphere_luhmann(
      double smooth_factor,
      double dipole_factor = 0.31, double csheet_factor = 0.15,
      double m_factor = earth_surface_avg_field);
  virtual ~msphere_luhmann() override;

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {return {0, 0, 0};};

  const metric_spherical* metric() const {return metric_;};
 private:
  const double c_bar_, d_bar_, idelta_;
  const metric_spherical* metric_;
};



} // end namespace gyronimo.

#endif // GYRONIMO_MSPHERE_LUHMANN
