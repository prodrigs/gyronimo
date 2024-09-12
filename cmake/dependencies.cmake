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
message(STATUS "  include: " ${GSL_INCLUDE_DIRS})
message(STATUS "  library: " ${GSL_LIBRARIES})

find_package(Boost 1.73.0 REQUIRED)
message(STATUS "  include: " ${Boost_INCLUDE_DIRS})

# if supporting VMEC, finds the needed resources (netcdf-cxx4 and netcdf):
if(SUPPORT_VMEC)
  message(STATUS "Configuring VMEC support (SUPPORT_VMEC=ON)")

  # tries to find the ncxx4-config utility in the current PATH:
  find_program(ncxx4-config_found "ncxx4-config")
  if(ncxx4-config_found)
    message(STATUS "  found ncxx4-config: " "${ncxx4-config_found}")
    execute_process(COMMAND "${ncxx4-config_found}" --includedir
        OUTPUT_VARIABLE NCXX4_INCLUDE_DIRS)
    string(STRIP ${NCXX4_INCLUDE_DIRS} NCXX4_INCLUDE_DIRS)
    execute_process(COMMAND "${ncxx4-config_found}" --libs
        OUTPUT_VARIABLE ncxx4_libs_raw)
    string(REGEX MATCH "^\ ?-L[^\ ]+" part1 ${ncxx4_libs_raw})  # removes...
    string(REGEX MATCH "\ -l[^\ ]+" part2 ${ncxx4_libs_raw})  # ...extraneous...
    string(JOIN " " NCXX4_LIBRARIES ${part1} ${part2})  # ... -lnetcdf info.
    message(STATUS "  include: " ${NCXX4_INCLUDE_DIRS})
    message(STATUS "  library: " ${NCXX4_LIBRARIES})
  else()
    message(FATAL_ERROR "ncxx4-config not found in PATH")
  endif()

  # tries to find the nc-config utility in the current PATH:
  find_program(nc-config_found "nc-config")
  if(nc-config_found)
    message(STATUS "  found nc-config: " "${nc-config_found}")
    execute_process(COMMAND "${nc-config_found}" --libs
        OUTPUT_VARIABLE nc_libs_raw)
    string(STRIP ${nc_libs_raw} NC_LIBRARIES)
    message(STATUS "  library: " ${NC_LIBRARIES})
  else()
    message(FATAL_ERROR "nc-config not found in PATH")
  endif()
else()
  message(STATUS "Skipping VMEC support (SUPPORT_VMEC=OFF)")
endif()
