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

# @dependencies.cmake, this file is part of ::gyronimo::

find_package(GSL REQUIRED)
  include_directories(${GSL_INCLUDE_DIRS})
  link_libraries(${GSL_LIBRARIES})
  message(STATUS "  include: " ${GSL_INCLUDE_DIRS})
  message(STATUS "  library: " ${GSL_LIBRARIES})

find_package(Boost 1.73.0 REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
  message(STATUS "  include: " ${Boost_INCLUDE_DIRS})

# add libraries to provide VMEC support (ncxx4 and dependencies) if required;
# otherwise, removes all related source files from the project:
if(SUPPORT_VMEC)
  message(STATUS "Configuring VMEC support (SUPPORT_VMEC=ON)")
  find_program(ncxx4_found "ncxx4-config")
  if(NOT ncxx4_found)
    message(FATAL_ERROR "ncxx4-config not found, install Unidata's netcdf-cxx4")
  endif()
  message(STATUS "  found ncxx4-config: " "${ncxx4_found}")
  execute_process(COMMAND "ncxx4-config" --includedir
      OUTPUT_VARIABLE ncxx4_include_dirs)
  string(STRIP ${ncxx4_include_dirs} ncxx4_include_dirs)
  execute_process(COMMAND "ncxx4-config" --libs OUTPUT_VARIABLE ncxx4_libs_raw)
  string(REGEX MATCH "^\ ?-L[^\ ]+" part1 ${ncxx4_libs_raw})
  string(REGEX MATCH "\ -l[^\ ]+" part2 ${ncxx4_libs_raw})
  string(JOIN " " ncxx4_libraries ${part1} ${part2})
  message(STATUS "  include: " ${ncxx4_include_dirs})
  message(STATUS "  library: " ${ncxx4_libraries})
  include_directories(${ncxx4_include_dirs})
  link_libraries(${ncxx4_libraries})
else()
  message(STATUS "Skipping VMEC support (SUPPORT_VMEC=OFF)")
  list(REMOVE_ITEM gyronimo_sources
      ${PROJECT_SOURCE_DIR}/gyronimo/metrics/metric_vmec.cc
      ${PROJECT_SOURCE_DIR}/gyronimo/parsers/parser_vmec.cc
      ${PROJECT_SOURCE_DIR}/gyronimo/metrics/morphism_vmec.cc
      ${PROJECT_SOURCE_DIR}/gyronimo/fields/equilibrium_vmec.cc)
  list(REMOVE_ITEM apps_sources
      ${PROJECT_SOURCE_DIR}/misc/apps/vmecdump.cc
      ${PROJECT_SOURCE_DIR}/misc/apps/vmectrace.cc)
endif()
