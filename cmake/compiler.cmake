# ::gyronimo:: - gyromotion for the people, by the people -
# An object-oriented library for gyromotion applications in plasma physics.
# Copyright (C) 2023 Paulo Rodrigues.

# ::gyronimo:: is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ::gyronimo:: is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

# @compiler.cmake, this file is part of ::gyronimo::

# detect platform:
string(COMPARE EQUAL "Darwin"  ${CMAKE_HOST_SYSTEM_NAME} OS_X)
string(COMPARE EQUAL "Linux"   ${CMAKE_HOST_SYSTEM_NAME} LINUX)
string(COMPARE EQUAL "Windows" ${CMAKE_HOST_SYSTEM_NAME} WINDOWS)
if(LINUX)
  set(OS_STRING "linux")
elseif(WINDOWS)
  set(OS_STRING "windows")
elseif(OS_X)
  set(OS_STRING "macOS")
else()
  set(OS_STRING "Unknown")
endif()

# compiler checkup and flags:
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# GNU section:
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  add_compile_options(-Wfatal-errors -Wpedantic)
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-Og -ggdb3)
  elseif(CMAKE_BUILD_TYPE MATCHES "Release")
    add_compile_options(-O2)  # overrides default -O3 optimisation flag
  elseif(CMAKE_BUILD_TYPE MATCHES "Plain")
# keep flag list clean here, allow it to be prefixed by -DCMAKE_CXX_FLAGS=...
  endif()
endif()
