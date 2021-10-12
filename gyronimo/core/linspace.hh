// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @linspace.hh

#ifndef GYRONIMO_LINSPACE
#define GYRONIMO_LINSPACE

#include <algorithm>
#include <gyronimo/core/error.hh>

namespace gyronimo {

//! Returns a `Container` of evenly spaced samples, similar to numpy's linspace.
template<typename Container> requires
  std::ranges::sized_range<Container> &&
  std::constructible_from<Container, size_t> &&
  std::floating_point<typename Container::value_type>
inline
Container linspace(
    const typename Container::value_type& start,
    const typename Container::value_type& end,
    size_t number) {
  if(end <= start || number < 2) error(
      __func__, __FILE__, __LINE__, "inconsistent arguments.", 1);
  Container samples(number);
  double delta = (end - start)/(number - 1);
  std::ranges::transform(
      std::views::iota(0u, number), std::ranges::begin(samples),
      [start, delta](size_t i) {return start + i*delta;});
  return samples;
}

} // end namespace gyronimo.

#endif // GYRONIMO_LINSPACE
