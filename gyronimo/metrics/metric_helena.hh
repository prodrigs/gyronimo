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

// @metric_helena.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_HELENA
#define GYRONIMO_METRIC_HELENA

#include <gyronimo/interpolators/interpolator2d.hh>
#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism_helena.hh>

namespace gyronimo {

//! Metric for `HELENA` field-aligned coordinates @f$\{s, \chi, \phi\}@f$.
/*!
    The contravariant coordinates are the square root of the poloidal flux per
    radian normalised to its value at the boundary (i.e.,
    @f$s=\sqrt{\Psi/\Psi_b}@f$), the angle on the poloidal cross section such
    that @f$B^\phi=q B^\chi@f$, with @f$q@f$ the safety factor, measured
    *counterclockwise* from the left-hand side midplane, and the toroidal angle
    measured *clockwise* when looking from the torus top. Both angles are in
    rads.

    The metric is built from the information provided by a `parser_helena`
    object inherited from the defining `morphism_helena`. The actual type of 2d
    interpolators to use is set by the specific `interpolator2d_factory` object
    pointer provided to the constructor. Notice that Christoffel symbols are
    computed using their `metric_covariant` implementation, which requires 1st
    order derivatives of the metric, rather than the `metric_connected` one to
    avoid resorting to the 2nd order derivatives of the `morphism_helena`
    interpolators.
*/
class metric_helena : public metric_connected {
 public:
  metric_helena(
      const morphism_helena* morph, const interpolator2d_factory* ifactory);
  virtual ~metric_helena() override;
  virtual SM3 operator()(const IR3& q) const override final;
  virtual dSM3 del(const IR3& q) const override final;
  virtual ddIR3 christoffel_first_kind(const IR3& q) const override final;
  virtual ddIR3 christoffel_second_kind(const IR3& q) const override final;

  const parser_helena* parser() const { return parser_; };
  const morphism_helena* my_morphism() const {
    return static_cast<const morphism_helena*>(metric_connected::my_morphism());
  };
 private:
  const parser_helena* parser_;
  interpolator2d *guu_, *guv_, *gvv_, *gww_;
  double R0_, squaredR0_;
};

inline ddIR3 metric_helena::christoffel_first_kind(const IR3& q) const {
  return this->gyronimo::metric_covariant::christoffel_first_kind(q);
}
inline ddIR3 metric_helena::christoffel_second_kind(const IR3& q) const {
  return this->gyronimo::metric_covariant::christoffel_second_kind(q);
}

}  // end namespace gyronimo

#endif  // GYRONIMO_METRIC_HELENA
