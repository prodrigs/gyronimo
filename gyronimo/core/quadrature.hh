// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @quadrature.hh

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
