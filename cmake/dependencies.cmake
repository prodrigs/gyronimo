# ::gyronimo:: - gyromotion for the people, by the people -
# An object-oriented library for gyromotion applications in plasma physics.
# Copyright (C) 2023-2024 Paulo Rodrigues.

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

# @dependencies.cmake, this file is part of ::gyronimo::

find_package(GSL REQUIRED)
message(STATUS "  include: " "${GSL_INCLUDE_DIRS}")
message(STATUS "  library: " "${GSL_LIBRARIES}")

find_package(Boost 1.73.0 REQUIRED)
message(STATUS "  include: " "${Boost_INCLUDE_DIRS}")

# add libraries to provide VMEC support (netcdf-cxx4) if required;
if(SUPPORT_VMEC)
  message(STATUS "Configuring VMEC support (SUPPORT_VMEC=ON)")
  find_path(ncxx4_include_dirs NAMES ncOpaqueType.h REQUIRED)
  find_library(ncxx4_libraries NAMES netcdf-cxx4 REQUIRED)
  message(STATUS "  include: " "${ncxx4_include_dirs}")
  message(STATUS "  library: " "${ncxx4_libraries}")
else()
  message(STATUS "Skipping VMEC support (SUPPORT_VMEC=OFF)")
endif()
