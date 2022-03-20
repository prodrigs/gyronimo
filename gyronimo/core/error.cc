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

// @error.cc, this file is part of ::gyronimo::

#include <cstdlib>
#include <iostream>
#include <gyronimo/core/error.hh>

namespace gyronimo {

//! Prints a warning message to `stderr`.
void warning(const char* message) {
  std::cerr << "gyronimo::warning: " << message << std::endl;
}

//! Prints an error message to `stderr` and exits execution with `exit_code`.
void error(
    const char* caller_name, const char* file_name,
    const int line_num, const char* message, const int exit_code) {
  std::cerr << "gyronimo::" << caller_name << ": " << message << std::endl;
  std::cerr << file_name << ", line " << line_num << ".\n";
  std::cerr << "error code: " << exit_code << ".\n";
  std::exit(exit_code);
}

} // end namespace gyronimo.
