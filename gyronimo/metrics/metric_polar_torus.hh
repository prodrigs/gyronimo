// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues and Manuel Assunção.

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

// @metric_polar_torus.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_POLAR_TORUS
#define GYRONIMO_METRIC_POLAR_TORUS

#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism_polar_torus.hh>

namespace gyronimo {

//! Metric for geometrical toroidal coordinates @f$\{r, \theta, \phi\}@f$.
/*!
    The contravariant coordinates are the shortest distance (normalized to the
    length `minor_radius`) to the circular line at `major_radius` from the torus
    axis of symmetry (which defines the torus' midplane), the angle measured
    counterclockwise on the poloidal cross section from the low-field side
    midplane, and the toroidal angle measured clockwise when looking from the
    torus' top. The lengths `minor_radius` and `major_radius` are in SI units
    and both angles are in rads.
*/
class metric_polar_torus : public metric_connected {
 public:
  metric_polar_torus(const morphism_polar_torus* morph);
  virtual ~metric_polar_torus() override {};
  virtual SM3 operator()(const IR3& q) const override final;
  virtual SM3 inverse(const IR3& q) const override final;
  virtual dSM3 del(const IR3& q) const override final;
  virtual double jacobian(const IR3& q) const override final;
  virtual IR3 to_covariant(const IR3& B, const IR3& q) const override final;
  virtual IR3 to_contravariant(const IR3& B, const IR3& q) const override final;

  double minor_radius() const { return minor_radius_; };
  double major_radius() const { return major_radius_; };
  double iaspect_ratio() const { return iaspect_ratio_; };
  const morphism_polar_torus* my_morphism() const {
    return static_cast<const morphism_polar_torus*>(
        metric_connected::my_morphism());
  };
 private:
  const double minor_radius_, major_radius_;
  const double minor_radius_squared_, iminor_radius_squared_;
  const double major_radius_squared_, iaspect_ratio_;
};

}  // end namespace gyronimo.

#endif  // GYRONIMO_METRIC_POLAR_TORUS
