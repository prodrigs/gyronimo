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

// @cubic_gsl.cc, this file is part of ::gyronimo::

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>

namespace gyronimo {

cubic_gsl::cubic_gsl(const dblock& x_range, const dblock& y_range)
    : spline1d_gsl() {
  spline_ = gsl_spline_alloc(gsl_interp_cspline, x_range.size());
  if (!spline_) error(
      __func__, __FILE__, __LINE__, "cannot allocate spline.", 1);
  gsl_spline_init(spline_, x_range.data(), y_range.data(), x_range.size());
}
cubic_gsl::~cubic_gsl() {
  if (spline_) delete spline_;
}

} // end namespace gyronimo.
