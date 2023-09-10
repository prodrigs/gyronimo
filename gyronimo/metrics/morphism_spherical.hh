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

// @morphism_spherical.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_SPHERICAL
#define GYRONIMO_MORPHISM_SPHERICAL

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Morphism from spherical coordinates @f$\{r, \phi, \theta\}@f$.
/*!
    The contravariant coordinates are the distance to the origin (normalised to
    `Lref` in SI units), the angle measured from the `z` axis (i.e., co-latitude
    measured from the north pole), and the angle measured from the `x` axis
    counterclockwise when seen from the north pole. Both angles are in rads.
*/
class morphism_spherical : public morphism {
 public:
  morphism_spherical(const double& Lref);
  virtual ~morphism_spherical() override {};

  virtual IR3 operator()(const IR3& q) const override final;
  virtual IR3 inverse(const IR3& x) const override final;
  virtual dIR3 del(const IR3& q) const override final;
  virtual ddIR3 ddel(const IR3& q) const override final;

  virtual double jacobian(const IR3& q) const override final;
  virtual dIR3 del_inverse(const IR3& q) const override final;

  double Lref() const { return Lref_; };
 private:
  const double Lref_, iLref_, Lref3_;
};

}  // end namespace gyronimo

#endif  // GYRONIMO_MORPHISM_SPHERICAL
