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

// @generators.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_GENERATORS
#define GYRONIMO_GENERATORS

#include <ranges>

namespace gyronimo {

template<typename Container>
concept SizedContiguousRange = 
  std::ranges::sized_range<Container> &&
  std::ranges::contiguous_range<Container>;

template<SizedContiguousRange Container> requires
  std::constructible_from<Container, size_t>
Container generate_sized(size_t size) {
  return Container(size);
}

template<SizedContiguousRange Container> requires
  (! std::constructible_from<Container, size_t>)
Container generate_sized(size_t size) {
  Container sample;
  if(size != sample.size()) gyronimo::error(__func__, __FILE__, __LINE__,
      "required size does not match template deduction.", 1);
  return sample;
}

} // end namespace gyronimo.

#endif // GYRONIMO_GENERATORS
