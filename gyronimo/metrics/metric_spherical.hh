// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues and Manuel Assunção.

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

//! Covariant metric for spherical coordinates.
/*!
    The three contravariant coordinates are the distance to the origin
    normalized to `Lref` (`u`, with `Lref` in SI), the polar angle
    measured from the z-axis (co-latitude `v`, in rads), and the azimuthal angle
    (`w`, also in rads) measured clockwise when looking from the origin along
    the z-axis. Some inherited methods are overriden for efficiency.
*/
class metric_spherical : public metric_connected {
 public:
  metric_spherical(const morphism_spherical* morph);
  virtual ~metric_spherical() override {};

  virtual SM3 operator()(const IR3& r) const override;
  virtual SM3 inverse(const IR3& r) const override;
  virtual dSM3 del(const IR3& r) const override;

  virtual double jacobian(const IR3& r) const override;
  virtual IR3 del_jacobian(const IR3& r) const override;
  virtual IR3 to_covariant(const IR3& B, const IR3& r) const override;
  virtual IR3 to_contravariant(const IR3& B, const IR3& r) const override;

  virtual ddIR3 christoffel_first_kind(const IR3& q) const override;
  virtual ddIR3 christoffel_second_kind(const IR3& q) const override;

  double Lref() const { return Lref_; };

  virtual const morphism_spherical* my_morphism() const override {
    return static_cast<const morphism_spherical*>(
        metric_connected::my_morphism());
  };
 private:
  const double Lref_, Lref_squared_;
  const double Lref_cube_, iLref_squared_;
};

}  // end namespace gyronimo.

#endif  // GYRONIMO_METRIC_SPHERICAL
