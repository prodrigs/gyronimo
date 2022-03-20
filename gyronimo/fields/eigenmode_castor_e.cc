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

// @eigenmode_castor_e.cc, this file is part of ::gyronimo::

#include <gyronimo/fields/eigenmode_castor_e.hh>

namespace gyronimo{

IR3 eigenmode_castor_e::covariant(const IR3& position, double time) const {
  IR3 A = A_.covariant(position, time);
  return {std::real(-iw_*A[IR3::u]),
      std::real(-iw_*A[IR3::v]), std::real(-iw_*A[IR3::w])};
}
IR3 eigenmode_castor_e::contravariant(const IR3& position, double time) const {
  IR3 A = A_.contravariant(position, time);
  return {std::real(-iw_*A[IR3::u]),
      std::real(-iw_*A[IR3::v]), std::real(-iw_*A[IR3::w])};
}

} // end namespace gyronimo.
