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

// @metric_cartesian.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_CARTESIAN
#define GYRONIMO_METRIC_CARTESIAN

#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism_cartesian.hh>

namespace gyronimo {

//! Trivial covariant metric for cartesian space.
class metric_cartesian : public metric_connected {
 public:
  metric_cartesian(const morphism_cartesian* m) : metric_connected(m) {};
  virtual ~metric_cartesian() override {};
  virtual SM3 operator()(const IR3& q) const override final;
  virtual SM3 inverse(const IR3& q) const override final;
  virtual dSM3 del(const IR3& q) const override final;
  virtual dSM3 del_inverse(const IR3& q) const override final;
  virtual double jacobian(const IR3& q) const override final;
  virtual IR3 del_jacobian(const IR3& q) const override final;
  virtual ddIR3 christoffel_first_kind(const IR3& q) const override final;
  virtual ddIR3 christoffel_second_kind(const IR3& q) const override final;
  virtual IR3 to_covariant(const IR3& B, const IR3& q) const override final;
  virtual IR3 to_contravariant(const IR3& B, const IR3& q) const override final;
  virtual IR3 inertial_force(
      const IR3& q, const IR3& dot_q) const override final;

  const morphism_cartesian* my_morphism() const {
    return static_cast<const morphism_cartesian*>(
        metric_connected::my_morphism());
  };
};

inline SM3 metric_cartesian::operator()(const IR3& q) const {
  return {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
}
inline SM3 metric_cartesian::inverse(const IR3& q) const {
  return {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
}
inline dSM3 metric_cartesian::del(const IR3& q) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
inline dSM3 metric_cartesian::del_inverse(const IR3& q) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
inline double metric_cartesian::jacobian(const IR3& q) const { return 1.0; }
inline IR3 metric_cartesian::del_jacobian(const IR3& q) const {
  return {0.0, 0.0, 0.0};
}
inline ddIR3 metric_cartesian::christoffel_first_kind(const IR3& q) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
inline ddIR3 metric_cartesian::christoffel_second_kind(const IR3& q) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
inline IR3 metric_cartesian::inertial_force(
    const IR3& q, const IR3& dot_q) const {
  return {0.0, 0.0, 0.0};
}
inline IR3 metric_cartesian::to_covariant(const IR3& B, const IR3& q) const {
  return B;
}
inline IR3 metric_cartesian::to_contravariant(
    const IR3& B, const IR3& q) const {
  return B;
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_METRIC_CARTESIAN
