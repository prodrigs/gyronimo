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

// @quadrature.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_QUADRATURE
#define GYRONIMO_QUADRATURE

#include <ranges>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptors.hpp>

namespace gyronimo {

//! Trapezoidal quadrature of a function uniformly sampled in [0:1].
template <typename Range> requires
  std::ranges::sized_range<Range> &&
  std::floating_point<typename Range::value_type> &&
  requires(const Range& c) {c.front(); c.back();}
Range::value_type trapezoidal(const Range& samples) {
  auto n = std::ranges::size(samples);
  auto boundary_correction = std::midpoint(samples.front(), samples.back());
  return boost::accumulate(samples, -boundary_correction)/(n - 1);
}

//! Trapezoidal quadrature of a function non-uniformly sampled.
template<typename RangeIn1, typename RangeIn2, typename RangeOut> requires
  std::floating_point<typename RangeIn1::value_type> &&
  std::floating_point<typename RangeIn2::value_type> &&
  std::floating_point<typename RangeOut::value_type> &&
  std::ranges::random_access_range<RangeIn1> &&
  std::ranges::random_access_range<RangeIn2> &&
  std::ranges::random_access_range<RangeOut> &&
  requires(const RangeOut& c) {c.back();}
RangeOut::value_type trapezoidal(
    const RangeIn1& samples, const RangeIn2& grid, RangeOut& quadrature) {
  auto q = std::ranges::begin(quadrature);
  auto y = std::ranges::begin(samples);
  auto x = std::ranges::begin(grid);
  q[0] = 0.0;
  for(size_t i = 1;i < std::ranges::size(samples); i++)
    q[i] = q[i - 1] + 0.5*(y[i] + y[i - 1])*(x[i] - x[i - 1]);
  return quadrature.back();
}

//! Simpson quadrature of a function uniformly sampled in [0:1].
template <typename Range> requires
  std::ranges::sized_range<Range> &&
  std::floating_point<typename Range::value_type> &&
  requires(const Range& c) {c.front(); c.back();}
Range::value_type simpson(const Range& samples) {
  using namespace boost::adaptors;
  auto n = std::ranges::size(samples);
  const typename Range::value_type zero = 0;
  auto sum_even = boost::accumulate(samples | strided(2), zero);
  auto sum_odd = boost::accumulate(samples | sliced(1, n) | strided(2), zero);
  auto boundary_correction = samples.front() + (n % 2 ?
      samples.back() : boost::inner_product(
          samples | sliced(n - 3, n),
          std::initializer_list<int>{1, -4, 11}, zero)/4);
  return (2*sum_even + 4*sum_odd - boundary_correction)/(3*(n - 1));
}

} // end namespace gyronimo.

#endif // GYRONIMO_QUADRATURE
