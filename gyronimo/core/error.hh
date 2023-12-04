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

// @error.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_ERROR
#define GYRONIMO_ERROR

#include <string>

namespace gyronimo {

void warning(const char* message);
void error(
    const char* caller_name, const char* file_name,
    const int line_num, const std::string& message, const int exit_code);

} // end namespace gyronimo.

#endif // GYRONIMO_ERROR.
