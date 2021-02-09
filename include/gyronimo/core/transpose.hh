// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @transpose.hh

#ifndef GYRONIMO_TRANSPOSE
#define GYRONIMO_TRANSPOSE

#include <ranges>
#include <gyronimo/core/error.hh>

namespace gyronimo {

//! Returns the transpose of a `nfast*nslow` matrix in a linear `Container`.
template<typename Container> requires 
  std::copy_constructible<Container> &&
  std::ranges::sized_range<Container> &&
  std::ranges::random_access_range<Container>
inline
Container transpose(const Container& original, size_t nfast) {
  if (std::size(original)%nfast) error(
        __func__, __FILE__, __LINE__, "inconsistent nfast and matrix size.", 1);
  auto nslow = std::size(original)/nfast;
  Container transposed(original);
  for (size_t p = 0; p < nslow*nfast; p++) {
    size_t i = p/nfast, j = p%nfast;
    transposed[i + j*nslow] = original[j + i*nfast];
  }
  return transposed;
}

} // end namespace gyronimo.

#endif // GYRONIMO_TRANSPOSE
