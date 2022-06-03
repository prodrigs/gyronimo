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

// @morphism_helena.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_HELENA
#define GYRONIMO_MORPHISM_HELENA

#include <gsl/gsl_vector.h>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/parsers/parser_helena.hh>
#include <gyronimo/interpolators/interpolator2d.hh>

namespace gyronimo {

//! Morphism from `HELENA` curvilinear coordinates.
class morphism_helena : public morphism {
 public:
  morphism_helena(
      const parser_helena *parser, const interpolator2d_factory *ifactory);
  virtual ~morphism_helena() override;
  virtual IR3 operator()(const IR3& q) const override;
  virtual IR3 inverse(const IR3& x) const override;
  virtual dIR3 del(const IR3& q) const override;
  const parser_helena* parser() const {return parser_;};
  double reduce_chi(double chi) const;
 private:
  static double reduce(double x, double l) {return (x -= l*std::floor(x/l));};
  const parser_helena *parser_;
  interpolator2d *R_, *z_;
  static int root_f(const gsl_vector *q, void *target, gsl_vector *f);
  struct root_target {double R, z;interpolator2d *i_R, *i_z;};
};

}

#endif // GYRONIMO_MORPHISM_HELENA
