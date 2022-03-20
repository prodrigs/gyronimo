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

// @bspline3_boost.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_BSPLINE3_BOOST
#define GYRONIMO_BSPLINE3_BOOST

#include <gyronimo/interpolators/interpolator1d.hh>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

namespace gyronimo {

//! Cubic B-spline by [boost](https://www.boost.org).
/*!
    If ommited from the constructor, `left_prime` and `right_prime` will be
    computed from the sampled data in `ordinates` using forward/backward finite
    differences (errors of the order of the grid size are implied, consider more
    accurate options).
*/
class bspline3_boost : public interpolator1d {
 public:
  bspline3_boost(const dblock& ordinates, double abcissa, double step)
    : spline_(ordinates.begin(), ordinates.end(), abcissa, step) {};
  bspline3_boost( const dblock& ordinates, double abcissa, double step,
      double left_prime, double right_prime)
    : spline_(ordinates.begin(), ordinates.end(),
          abcissa, step, left_prime, right_prime) {};
  virtual ~bspline3_boost() final {};
  double operator()(double x) const final {return spline_(x);};
  double derivative(double x) const final {return spline_.prime(x);};
  double derivative2(double x) const final {return spline_.double_prime(x);};
 private:
  boost::math::interpolators::cardinal_cubic_b_spline<double> spline_;
};

//! Factory for cubic B-spline by [boost](https://www.boost.org).
/*!
    Builds interpolators according to a specified policy to get the left and
    right derivatives: native uses the boost implementation (forward/backward
    finite differences), periodic assumes `ordinates.front() ==
    ordinates.back()` and uses finite differences centred at the first and last
    samples, `periodic_short` assumes periodic repetition of the sampled data
    (i.e., `ordinates.front() == (*ordinates.end()`).
*/
class bspline3_boost_factory : public interpolator1d_factory {
 public:
  enum policy {native, periodic, periodic_short};
  bspline3_boost_factory(const policy p) : policy_(p) {};
  virtual interpolator1d* interpolate_data(
      const dblock& abcissas, const dblock& ordinates) const final;
 private:
  const policy policy_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_BSPLINE3_BOOST
