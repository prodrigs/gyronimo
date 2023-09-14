// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Paulo Rodrigues and Manuel Assunção.

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

// @metric_spherical.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_SPHERICAL
#define GYRONIMO_METRIC_SPHERICAL

#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism_spherical.hh>

namespace gyronimo {

//! Metric for spherical coordinates @f$\{r, \phi, \theta\}@f$.
/*!
    The contravariant coordinates are the distance to the origin (normalised to
    `Lref` in SI units), the angle measured from the `z` axis (i.e., co-latitude
    measured from the north pole), and the angle measured from the `x` axis
    counterclockwise when seen from the north pole. Both angles are in rads.
*/
class metric_spherical : public metric_connected {
 public:
  metric_spherical(const morphism_spherical* morph);
  virtual ~metric_spherical() override {};
  virtual SM3 operator()(const IR3& q) const override;
  virtual SM3 inverse(const IR3& q) const override;
  virtual dSM3 del(const IR3& q) const override;
  virtual double jacobian(const IR3& q) const override;
  virtual IR3 del_jacobian(const IR3& q) const override;
  virtual IR3 to_covariant(const IR3& B, const IR3& q) const override;
  virtual IR3 to_contravariant(const IR3& B, const IR3& q) const override;
  virtual ddIR3 christoffel_first_kind(const IR3& q) const override;
  virtual ddIR3 christoffel_second_kind(const IR3& q) const override;
  virtual IR3 inertial_force(const IR3& q, const IR3& dot_q) const override;

  double Lref() const { return Lref_; };
  const morphism_spherical* my_morphism() const {
    return static_cast<const morphism_spherical*>(
        metric_connected::my_morphism());
  };
 private:
  const double Lref_, Lref2_, Lref3_, iLref2_;
};

}  // end namespace gyronimo.

#endif  // GYRONIMO_METRIC_SPHERICAL
