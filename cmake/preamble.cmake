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

# @preamble.cmake, this file is part of ::gyronimo::

# disallows in-source build:
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR “In-source build not allowed.”)
endif()

# checks and sets build type (defaults to Release):
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
string(FIND "Debug Release Plain" ${CMAKE_BUILD_TYPE} KNOWN_BUILD_TYPE)
if(${KNOWN_BUILD_TYPE} LESS "0")
  message(FATAL_ERROR "Unknown build type " ${CMAKE_BUILD_TYPE} ".")
endif()

# gets the latest commit hash and routes it to version.hh via compile flag:
find_program(git_found "git")
if(git_found)
  execute_process(
    COMMAND git rev-parse --short=16 HEAD WORKING_DIRECTORY
    ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(GIT_COMMIT_HASH "unavailable")
endif()
add_compile_definitions("-DGYRONIMO_GIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")

# extracts version numbers from version.hh, stores output in gyronimo_version:
set(version_major_regex "[ \t]*constexpr int version_major = ([0-9]+);")
set(version_minor_regex "[ \t]*constexpr int version_minor = ([0-9]+);")
set(version_patch_regex "[ \t]*constexpr int version_patch = ([0-9]+);")
file(STRINGS ${CMAKE_SOURCE_DIR}/gyronimo/version.hh
  string_to_replace REGEX ${version_major_regex})
string(REGEX REPLACE
  ${version_major_regex} "\\1" version_major ${string_to_replace})
file(STRINGS ${CMAKE_SOURCE_DIR}/gyronimo/version.hh
  string_to_replace REGEX ${version_minor_regex})
string(REGEX REPLACE
  ${version_minor_regex} "\\1" version_minor ${string_to_replace})
file(STRINGS ${CMAKE_SOURCE_DIR}/gyronimo/version.hh
  string_to_replace REGEX ${version_patch_regex})
string(REGEX REPLACE
  ${version_patch_regex} "\\1" version_patch ${string_to_replace})

string(APPEND gyronimo_version
  ${version_major} "." ${version_minor} "." ${version_patch})
