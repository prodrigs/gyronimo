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

// @codata.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CODATA
#define GYRONIMO_CODATA

namespace gyronimo {

//! Physical constants (SI units) by [CODATA](https://physics.nist.gov/cuu/Constants/index.html).
namespace codata {

constexpr double c = 2.99792458e8;
constexpr double e = 1.602176634e-19;
constexpr double mu0 = 1.25663706212e-6;
constexpr double m_alpha = 6.6446573357e-27;
constexpr double m_proton = 1.67262192369e-27;
constexpr double m_neutron = 1.67492749804e-27;
constexpr double m_electron = 9.1093837015e-31;

} // end namespace codata

} // end namespace gyronimo.

#endif // GYRONIMO_CODATA
