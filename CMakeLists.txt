# ::gyronimo:: - gyromotion for the people, by the people -
# An object-oriented library for gyromotion applications in plasma physics.
# Copyright (C) 2021-2022 Jorge Ferreira and Paulo Rodrigues.

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

cmake_minimum_required(VERSION 3.16)

# disallows in-source build:
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
message(FATAL_ERROR “In-source build not allowed.”)
endif()

# extracts the version numbers from version.hh:
set(version_major_regex "[ \t]*constexpr int version_major = ([0-9]+);")
set(version_minor_regex "[ \t]*constexpr int version_minor = ([0-9]+);")
file(STRINGS ${CMAKE_SOURCE_DIR}/gyronimo/version.hh
  string_to_replace REGEX ${version_major_regex})
string(REGEX REPLACE
  ${version_major_regex} "\\1" version_major ${string_to_replace})
file(STRINGS ${CMAKE_SOURCE_DIR}/gyronimo/version.hh
  string_to_replace REGEX ${version_minor_regex})
string(REGEX REPLACE
  ${version_minor_regex} "\\1" version_minor ${string_to_replace})
string(APPEND version_string ${version_major} "." ${version_minor})

# defaults to Release build type:
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message("*  Configuring gyronimo-"
  ${version_string} " in " ${CMAKE_BUILD_TYPE} " mode:")

# build options
option(BUILD_VMEC_SUPPORT "Add support to VMEC equilibria" OFF)

# detect platform
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

# sets project properties:
project(
  gyronimo
  VERSION ${version_string}
  DESCRIPTION
    "An object-oriented library for gyromotion applications in plasma physics."
  HOMEPAGE_URL https://github.com/prodrigs/gyronimo
  LANGUAGES CXX
)
include(GNUInstallDirs)
include(ExternalProject)

# compiler checkup and flags:
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  add_compile_options(-Wfatal-errors -Wpedantic)
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-Og -ggdb3)
  endif()
  if(CMAKE_BUILD_TYPE MATCHES "Release")
    add_compile_options(-march=native)
  endif()
endif()

# find all gyronimo sources
file(GLOB gyronimo_sources ${PROJECT_SOURCE_DIR}/gyronimo/*/*.cc)
# find all apps sources
file(GLOB apps_sources ${PROJECT_SOURCE_DIR}/misc/apps/*.cc)

# required libraries
find_package(GSL REQUIRED)
  include_directories(${GSL_INCLUDE_DIRS})
  set(gyronimo_depending_libs ${gyronimo_depending_libs} ${GSL_LIBRARIES})
find_package(Boost REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
# add libraries needed to provide VMEC support
if(BUILD_VMEC_SUPPORT)
  # find_package(PkgConfig REQUIRED)
  # pkg_check_modules(netcdf-cxx4 REQUIRED IMPORTED_TARGET netcdf-cxx4)
  find_program(NCXX4_CONFIG "ncxx4-config")
  if (NCXX4_CONFIG)
    execute_process(COMMAND "ncxx4-config" --libs
        OUTPUT_VARIABLE netcdf_cxx_libs)
    execute_process(COMMAND "ncxx4-config" --cflags
        OUTPUT_VARIABLE netcdf_cxx_flags)
  else()
    message(SEND_ERROR "ncxx4-config not found. Please install Unidata's netCDF-CXX4 package")
  endif()
  string(STRIP ${netcdf_cxx_libs} netcdf_cxx_libs)
  string(STRIP ${netcdf_cxx_flags} netcdf_cxx_flags)
else()
    # remove all VMEC related gyronimo and apps sources
    message(STATUS "Not building VMEC support!")
    list(REMOVE_ITEM gyronimo_sources
      ${PROJECT_SOURCE_DIR}/gyronimo/metrics/metric_vmec.cc
      ${PROJECT_SOURCE_DIR}/gyronimo/parsers/parser_vmec.cc
    )
    list(REMOVE_ITEM apps_sources
      ${PROJECT_SOURCE_DIR}/misc/apps/vmecdump.cc
      ${PROJECT_SOURCE_DIR}/misc/apps/vmectrace.cc
    )
endif()

# builds gyronimo library:
add_library(gyronimo SHARED ${gyronimo_sources})
set_target_properties(gyronimo PROPERTIES VERSION ${PROJECT_VERSION})
if(NOT BUILD_VMEC_SUPPORT)
  target_include_directories(gyronimo PUBLIC ${PROJECT_SOURCE_DIR})
  target_link_libraries(gyronimo PUBLIC ${gyronimo_depending_libs})
else()
  # target_link_libraries(gyronimo PUBLIC ${gyronimo_depending_libs} PkgConfig::netcdf-cxx4)
  target_include_directories(gyronimo PUBLIC ${PROJECT_SOURCE_DIR} ${netcdf_cxx_flags})
  target_link_libraries(gyronimo PUBLIC ${gyronimo_depending_libs} ${netcdf_cxx_libs})
endif()

# builds the apps:
foreach(target_file IN LISTS apps_sources)
  get_filename_component(target ${target_file} NAME_WE)
  add_executable(${target} ${target_file})
  target_include_directories(${target}
    PUBLIC ${PROJECT_SOURCE_DIR}/misc/apps/include)
  target_link_libraries(${target} PRIVATE gyronimo)
  install(TARGETS ${target} DESTINATION bin)
endforeach()

# builds documentation with doxygen (requires explicit `--target doc`):
find_package(Doxygen)
if(DOXYGEN_FOUND)
# set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "readme.md")
  set(DOXYGEN_PROJECT_NAME "::gyronimo::")
  set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
  set(DOXYGEN_GENERATE_TREEVIEW NO)
  set(DOXYGEN_DOT_TRANSPARENT YES)
  set(DOXYGEN_SHOW_USED_FILES NO)
  set(DOXYGEN_SHOW_NAMESPACES NO)
  set(DOXYGEN_DISABLE_INDEX YES)
  set(DOXYGEN_SEARCHENGINE NO)
  set(DOXYGEN_SHOW_FILES YES)
  set(DOXYGEN_RECURSIVE YES)
  set(DOXYGEN_WARNINGS YES)
  set(DOXYGEN_WARN_IF_DOC_ERROR YES)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
  set(DOXYGEN_QUIET YES)
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_MATHJAX_FORMAT SVG)
  set(DOXYGEN_MATHJAX_RELPATH "https://cdn.jsdelivr.net/npm/mathjax@2")
  set(DOXYGEN_HTML_OUTPUT doc-html)
  doxygen_add_docs(doc
    ${PROJECT_SOURCE_DIR}/gyronimo
    COMMENT "Extracting documentation with doxygen")
  install(
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc-html/
    TYPE DOC MESSAGE_NEVER
    OPTIONAL)
endif()

# installs gyronimo library:
install(TARGETS gyronimo LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/gyronimo
  TYPE INCLUDE MESSAGE_NEVER
  FILES_MATCHING PATTERN "*.hh")

# pkg-config support
#configure_file(
#  ${PROJECT_SOURCE_DIR}/misc/templates/gyronimo.pc.in gyronimo.pc @ONLY)
#install(FILES ${CMAKE_BINARY_DIR}/gyronimo.pc
#  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
