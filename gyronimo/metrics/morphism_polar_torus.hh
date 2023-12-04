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

// @morphism_polar_torus.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_POLAR_TORUS
#define GYRONIMO_MORPHISM_POLAR_TORUS

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Morphism from geometrical toroidal coordinates @f$\{r, \theta, \phi\}@f$.
/*!
    The contravariant coordinates are the shortest distance (normalized to the
    length `minor_radius`) to the circular line at `major_radius` from the torus
    axis of symmetry (which defines the torus' midplane), the angle measured
    counterclockwise on the poloidal cross section from the low-field side
    midplane, and the toroidal angle measured clockwise when looking from the
    torus' top. The lengths `minor_radius` and `major_radius` are in SI units
    and both angles are in rads.
*/
class morphism_polar_torus : public morphism {
 public:
  morphism_polar_torus(const double minor_radius, double major_radius);
  virtual ~morphism_polar_torus() override {};

  virtual IR3 operator()(const IR3& q) const override;
  virtual IR3 inverse(const IR3& x) const override;
  virtual dIR3 del(const IR3& q) const override;
  virtual ddIR3 ddel(const IR3& q) const override;

  virtual double jacobian(const IR3& q) const override;
  virtual dIR3 del_inverse(const IR3& q) const override;

  double minor_radius() const { return minor_radius_; };
  double major_radius() const { return major_radius_; };
  double iaspect_ratio() const { return iaspect_ratio_; };
 private:
  const double minor_radius_, major_radius_;
  const double iaspect_ratio_, volume_factor_;
  const double iminor_radius_;
};

}  // end namespace gyronimo

#endif  // GYRONIMO_MORPHISM_POLAR_TORUS
