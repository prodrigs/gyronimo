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

// @spline1d_gsl.cc, this file is part of ::gyronimo::

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/spline1d_gsl.hh>

namespace gyronimo {

spline1d_gsl::spline1d_gsl() : spline_(nullptr), acc_(nullptr) {
  acc_ = gsl_interp_accel_alloc();
}
spline1d_gsl::~spline1d_gsl() {
  if (acc_) delete(acc_);
}
double spline1d_gsl::operator()(double x) const {
  return gsl_spline_eval(spline_, x, acc_);
}
double spline1d_gsl::derivative(double x) const {
  return  gsl_spline_eval_deriv(spline_, x, acc_);
}
double spline1d_gsl::derivative2(double x) const {
  return  gsl_spline_eval_deriv2(spline_, x, acc_);
}

} // end namespace gyronimo.
