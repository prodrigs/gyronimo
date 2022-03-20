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

// @spline1d_gsl.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_SPLINE1D_GSL
#define GYRONIMO_SPLINE1D_GSL

#include <gsl/gsl_spline.h>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo {

//! Base class for 1d splines by [GSL](https://www.gnu.org/software/gsl).
class spline1d_gsl : public interpolator1d {
 public:
  spline1d_gsl();
  virtual ~spline1d_gsl() override;

  double operator()(double x) const final;
  double derivative(double x) const final;
  double derivative2(double x) const final;

 protected:
  gsl_spline *spline_;
  gsl_interp_accel *acc_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_SPLINE1D_GSL
