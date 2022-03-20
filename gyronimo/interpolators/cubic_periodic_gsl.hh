// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @cubic_periodic_gsl.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CUBIC_PERIODIC_GSL
#define GYRONIMO_CUBIC_PERIODIC_GSL

#include <gyronimo/interpolators/spline1d_gsl.hh>

namespace gyronimo {

//! Periodic cubic spline by [GSL](https://www.gnu.org/software/gsl).
class cubic_periodic_gsl : public spline1d_gsl {
 public:
  cubic_periodic_gsl(const dblock& x_range, const dblock& y_range);
  virtual ~cubic_periodic_gsl() final;
};

class cubic_periodic_gsl_factory : public interpolator1d_factory {
 public:
  virtual interpolator1d* interpolate_data(
      const dblock& x_range, const dblock& y_range) const final {
    return new cubic_periodic_gsl(x_range, y_range);
  };
};

} // end namespace gyronimo.

#endif // GYRONIMO_CUBIC_PERIODIC_GSL
