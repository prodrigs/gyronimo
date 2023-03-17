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

// @morphism_helena.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_HELENA
#define GYRONIMO_MORPHISM_HELENA

#include <gyronimo/interpolators/interpolator2d.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/parsers/parser_helena.hh>

namespace gyronimo {

//! Morphism for `HELENA` curvilinear coordinates.
/*!
    Builds the morphism from the information provided by a `parser_helena`
    object. The actual type of 2d interpolators to use is set by the specific
    `interpolator2d_factory` object pointer provided to the constructor.

    The right-handed curvilinear coordinates are the square root of the poloidal
    flux per radian normalised to its value at the boundary (`u`, or `HELENA`
    @f$s=\sqrt{\Psi/\Psi_b}@f$), an angle on the poloidal cross section (`v`, or
    `HELENA` @f$\chi : B^\phi=q B^\chi@f$, in rads) measured *counterclockwise*
    from the left-hand side midplane (the `x` axis), and the toroidal angle
    (`w`, or `HELENA` @f$\phi@f$, in rads) measured *clockwise* when looking
    from the torus top.
*/
class morphism_helena : public morphism {
 public:
  morphism_helena(
      const parser_helena* parser, const interpolator2d_factory* ifactory);
  virtual ~morphism_helena() override;
  virtual IR3 operator()(const IR3& q) const override;
  virtual IR3 inverse(const IR3& x) const override;
  virtual dIR3 del(const IR3& q) const override;
  virtual ddIR3 ddel(const IR3& q) const override;

  virtual double jacobian(const IR3& q) const override;
  virtual IR3 translation(const IR3& q, const IR3& delta) const override;
  const parser_helena* parser() const { return parser_; };
 private:
  const parser_helena* parser_;
  interpolator2d *R_, *z_;
  std::pair<double, double> reflection_past_axis(double s, double chi) const;
};

}  // namespace gyronimo

#endif  // GYRONIMO_MORPHISM_HELENA
