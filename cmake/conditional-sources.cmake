# ::gyronimo:: - gyromotion for the people, by the people -
# An object-oriented library for gyromotion applications in plasma physics.
# Copyright (C) 2024 Paulo Rodrigues.

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

# @conditional-sources.cmake, this file is part of ::gyronimo::

if(NOT SUPPORT_VMEC)
  list(REMOVE_ITEM gyronimo_sources
      ${PROJECT_SOURCE_DIR}/gyronimo/metrics/metric_vmec.cc
      ${PROJECT_SOURCE_DIR}/gyronimo/parsers/parser_vmec.cc
      ${PROJECT_SOURCE_DIR}/gyronimo/metrics/morphism_vmec.cc
      ${PROJECT_SOURCE_DIR}/gyronimo/fields/equilibrium_vmec.cc)
  list(REMOVE_ITEM apps_sources
      ${PROJECT_SOURCE_DIR}/misc/apps/vmecdump.cc
      ${PROJECT_SOURCE_DIR}/misc/apps/vmectrace.cc)
endif()
