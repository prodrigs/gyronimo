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

// @io_streaming.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_IO_STREAMING
#define GYRONIMO_IO_STREAMING

#include <iostream>

namespace gyronimo {

//! Read helper for indexed objects.
template<class T>
std::istream& operator>>(std::istream& is,T& x) {
  for (size_t i = 0; i < x.size(); i++) is >> x[i];
  return is;
}

//! Write helper for indexed objects.
template<class T>
std::ostream& operator<<(std::ostream& os,T& x) {
  for (size_t i = 0; i < x.size(); i++) os << x[i];
  return os;
}

} // end namespace gyronimo.

#endif // GYRONIMO_IO_STREAMING
