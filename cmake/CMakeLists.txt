# ::gyronimo:: - gyromotion for the people, by the people -
# An object-oriented library for gyromotion applications in plasma physics.
# Copyright (C) 2021-2023 Jorge Ferreira and Paulo Rodrigues.

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

# @CMakeLists.txt, this file is part of ::gyronimo::

cmake_minimum_required(VERSION 3.19)
include(${CMAKE_SOURCE_DIR}/cmake/preamble.cmake)
message("*  Configuring ::gyronimo::v"
  ${gyronimo_version} " in " ${CMAKE_BUILD_TYPE} " mode:")
project(
  gyronimo
  VERSION ${gyronimo_version}
  DESCRIPTION
    "An object-oriented library for gyromotion applications in plasma physics."
  HOMEPAGE_URL https://github.com/prodrigs/gyronimo
  LANGUAGES CXX
)
include(GNUInstallDirs)
include(ExternalProject)
include(${CMAKE_SOURCE_DIR}/cmake/compiler.cmake)

# define build options and their default values:
option(SUPPORT_VMEC "Add support to VMEC equilibria" OFF)

# check required dependencies, updates include and lib directories:
include(${CMAKE_SOURCE_DIR}/cmake/dependencies.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/doxygen.cmake)

# defines source/build/install for gyronimo library:
file(GLOB gyronimo_sources ${PROJECT_SOURCE_DIR}/gyronimo/*/*.cc)
add_library(gyronimo SHARED ${gyronimo_sources})
set_target_properties(gyronimo PROPERTIES VERSION ${PROJECT_VERSION})
install(TARGETS gyronimo LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/gyronimo
  TYPE INCLUDE MESSAGE_NEVER FILES_MATCHING PATTERN "*.hh")

# defines source/build/install for the available apps:
file(GLOB apps_sources ${PROJECT_SOURCE_DIR}/misc/apps/*.cc)
foreach(target_file IN LISTS apps_sources)
  get_filename_component(target ${target_file} NAME_WE)
  add_executable(${target} ${target_file})
  target_include_directories(${target}
      PUBLIC ${PROJECT_SOURCE_DIR}/misc/apps/include)
  target_link_libraries(${target} PRIVATE gyronimo)
  install(TARGETS ${target} DESTINATION bin)
endforeach()
