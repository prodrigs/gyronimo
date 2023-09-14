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

// @metric_cylindrical.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_CYLINDRICAL
#define GYRONIMO_METRIC_CYLINDRICAL

#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism_cylindrical.hh>

namespace gyronimo {

//! Metric for cylindrical coordinates @f$\{r, \phi, z\}@f$.
/*!
    The contravariant coordinates are the distance to the axis (normalised to
    `Lref` in SI units), the angle measured counterclockwise when seen from the
    top of the axis (in rads), and the length measured along the latter (also
    normalised to `Lref`).
*/
class metric_cylindrical : public metric_connected {
 public:
  metric_cylindrical(const morphism_cylindrical* morph);
  virtual ~metric_cylindrical() override {};
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

  double Lref() { return Lref_; };
  const morphism_cylindrical* my_morphism() const {
    return static_cast<const morphism_cylindrical*>(
        metric_connected::my_morphism());
  };
 private:
  const double Lref_, Lref2_, iLref2_, Lref3_;
};

inline SM3 metric_cylindrical::operator()(const IR3& q) const {
  return {Lref2_, 0.0, 0.0, Lref2_ * q[IR3::u] * q[IR3::u], 0.0, Lref2_};
}
inline SM3 metric_cylindrical::inverse(const IR3& q) const {
  return {
      iLref2_, 0.0, 0.0, iLref2_ / (q[IR3::u] * q[IR3::u]), 0.0, iLref2_};
}
inline double metric_cylindrical::jacobian(const IR3& q) const {
  return Lref3_ * q[IR3::u];
}
inline IR3 metric_cylindrical::del_jacobian(const IR3& q) const {
  return {Lref3_, 0.0, 0.0};
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_METRIC_CYLINDRICAL
