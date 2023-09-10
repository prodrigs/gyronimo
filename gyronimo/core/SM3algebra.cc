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

// @SM3algebra.cc, this file is part of ::gyronimo::

#include <gyronimo/core/SM3algebra.hh>

namespace gyronimo {

SM3 inverse(const SM3& m) {
  double ideterminant = 1.0/(
      m[SM3::uw]*m[SM3::uw]*m[SM3::vv] - 2.0*m[SM3::uv]*m[SM3::uw]*m[SM3::vw] +
      m[SM3::uv]*m[SM3::uv]*m[SM3::ww] + m[SM3::uu]*(
          m[SM3::vw]*m[SM3::vw] - m[SM3::vv]*m[SM3::ww]));
  return {
    ideterminant*(m[SM3::vw]*m[SM3::vw] - m[SM3::vv]*m[SM3::ww]),
    ideterminant*(m[SM3::uv]*m[SM3::ww] - m[SM3::uw]*m[SM3::vw]),
    ideterminant*(m[SM3::uw]*m[SM3::vv] - m[SM3::uv]*m[SM3::vw]),
    ideterminant*(m[SM3::uw]*m[SM3::uw] - m[SM3::uu]*m[SM3::ww]),
    ideterminant*(m[SM3::uu]*m[SM3::vw] - m[SM3::uv]*m[SM3::uw]),
    ideterminant*(m[SM3::uv]*m[SM3::uv] - m[SM3::uu]*m[SM3::vv])};
}

} // end namespace gyronimo.
