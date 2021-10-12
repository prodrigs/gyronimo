# ::gyronimo:: - gyromotion for the people, by the people -
# An object-oriented library for gyromotion applications in plasma physics.
# Copyright (C) 2021 Jorge Ferreira, Paulo Rodrigues.

# @CMakeLists.txt

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

# required libraries
find_package(GSL REQUIRED)
  include_directories(${GSL_INCLUDE_DIRS})
  set(gyronimo_depending_libs ${gyronimo_depending_libs} ${GSL_LIBRARIES})
find_package(Boost REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})

# builds gyronimo library:
file(GLOB gyronimo_sources ${PROJECT_SOURCE_DIR}/gyronimo/*/*.cc)
add_library(gyronimo SHARED ${gyronimo_sources})
set_target_properties(gyronimo PROPERTIES VERSION ${PROJECT_VERSION})
target_include_directories(gyronimo PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(gyronimo PUBLIC ${gyronimo_depending_libs})

# builds the apps:
file(GLOB apps_sources ${PROJECT_SOURCE_DIR}/misc/apps/*.cc)
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
