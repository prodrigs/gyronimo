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

// @interpolator2d.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_INTERPOLATOR2D
#define GYRONIMO_INTERPOLATOR2D

#include <gyronimo/core/dblock.hh>

namespace gyronimo {

//! Access interface for 2d interpolators.
/*!
    Notice that this class only requires the **access** functionality to be
    implemented (evaluation, derivatives, second derivatives). Check the
    documentation of `interpolator1d` and `interpolator1d_factory` for details
    about the creation of specific interpolator objects by abstract code.
*/
class interpolator2d {
 public:
  virtual ~interpolator2d() {};
  virtual double operator()(double x, double y) const = 0;
  virtual double partial_u(double x, double y) const = 0;
  virtual double partial_v(double x, double y) const = 0;
  virtual double partial2_uu(double x, double y) const = 0;
  virtual double partial2_uv(double x, double y) const = 0;
  virtual double partial2_vv(double x, double y) const = 0;
};

//! Creation interface for 2d interpolators.
/*!
    Check the documentation of interpolator1d and interpolator1d_factory for
    details about the creation of specific interpolator objects by abstract
    code.
*/
class interpolator2d_factory {
 public:
  virtual interpolator2d* interpolate_data(
      const dblock& x_range,
      const dblock& y_range,
      const dblock& z_range) const = 0;
};

} // end namespace gyronimo.

#endif // GYRONIMO_INTERPOLATOR2D
