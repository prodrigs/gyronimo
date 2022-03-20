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

// @transpose.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_TRANSPOSE
#define GYRONIMO_TRANSPOSE

#include <ranges>
#include <vector>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/dblock.hh>

namespace gyronimo {

//! Returns the transpose of a `nfast*nslow` matrix in a linear `Range`.
//@todo relace copy_constructible by something like constructible_by_size in
// order to avoid construction without self allocation.
template<typename Range> requires 
  std::copy_constructible<Range> &&
  std::ranges::sized_range<Range> &&
  std::ranges::random_access_range<Range>
inline
Range transpose(const Range& original, size_t nfast) {
  if (std::size(original)%nfast) error(
        __func__, __FILE__, __LINE__, "inconsistent nfast and matrix size.", 1);
  auto nslow = std::size(original)/nfast;
  Range transposed(original);
  for (size_t p = 0; p < nslow*nfast; p++) {
    size_t i = p/nfast, j = p%nfast;
    transposed[i + j*nslow] = original[j + i*nfast];
  }
  return transposed;
}

//! Returns the transpose of a `nfast*nslow` matrix wrapped inside a dblock.
//@todo remove the dependency on std::vector and ensure proper ownership.
//@todo centralise the actual transposing code to avoid redundant code.
inline
dblock_adapter<std::vector<double>> transpose(
    const dblock& original, size_t nfast) {
  auto nslow = std::size(original)/nfast;
  if (std::size(original)%nfast) error(
        __func__, __FILE__, __LINE__, "inconsistent nfast and matrix size.", 1);
  std::vector<double> transposed(original.begin(), original.end());
  for (size_t p = 0; p < nslow*nfast; p++) {
    size_t i = p/nfast, j = p%nfast;
    transposed[i + j*nslow] = original.data()[j + i*nfast];
  }
  return dblock_adapter<std::vector<double>>(std::move(transposed));
}

} // end namespace gyronimo.

#endif // GYRONIMO_TRANSPOSE
