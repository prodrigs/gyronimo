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

// @bspline3_boost.cc, this file is part of ::gyronimo::

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/bspline3_boost.hh>

namespace gyronimo {

interpolator1d* bspline3_boost_factory::interpolate_data(
    const dblock& x_range, const dblock& ordinates) const {
  if (x_range.size() != ordinates.size())
      error(__func__, __FILE__, __LINE__, "size mismatch", 1);
  double init_x = x_range.front();
  double step = (x_range.back() - x_range.front())/(x_range.size() - 1);
  if (policy_ == native)
    return new bspline3_boost(ordinates, init_x, step);
  else {
    size_t n = ordinates.size();
    double inv_double_step = 0.5/step;
    if (policy_ == periodic) {
      double left_right_prime =
          (ordinates[1] - ordinates[n - 2])*inv_double_step;
      return new bspline3_boost(
          ordinates, init_x, step, left_right_prime, left_right_prime);
    } else {
        if (policy_ == periodic_short) {
        double left_prime = (ordinates[1] - ordinates[n - 1])*inv_double_step;
        double right_prime = (ordinates[0] - ordinates[n - 2])*inv_double_step;
        return new bspline3_boost(
            ordinates, init_x, step, left_prime, right_prime);
      }
    }
  }
  error(__func__, __FILE__, __LINE__, "unknown policy", 1);
  return nullptr;
}

} // end namespace gyronimo.
