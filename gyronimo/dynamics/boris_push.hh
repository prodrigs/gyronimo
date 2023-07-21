// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Manuel Assunção and Paulo Rodrigues.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @boris_push.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_BORIS_PUSH
#define GYRONIMO_BORIS_PUSH

#include <gyronimo/core/IR3algebra.hh>

namespace gyronimo {

//! Boris cartesian-velocity pusher (magnetic field only).
IR3 boris_push(
    const IR3& velocity, const double& tildeOref, const double& B, const IR3& b,
    const double& dt);

//! Boris velocity pusher (electric and magnetic fields).
IR3 boris_push(
    const IR3& velocity, const double& tildeOref, const double& tildeEref,
    const IR3& E, const double& B, const IR3& b, const double& dt);

}  // end namespace gyronimo

#endif  // GYRONIMO_BORIS_PUSH
