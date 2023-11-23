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

# @doxygen.cmake, this file is part of ::gyronimo::

# builds documentation with doxygen (requires explicit `--target doc`):
if(SUPPORT_DOXYGEN)
  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Cannot configure doxygen (explicit -DSUPPORT_DOXYGEN=ON)")
  endif()
  message(STATUS "Configuring doxygen")
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
