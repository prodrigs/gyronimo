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

// @boris_push.cc, this file is part of ::gyronimo::

#include <gyronimo/core/contraction.hh>
#include <gyronimo/dynamics/boris_push.hh>

#include <cmath>

namespace gyronimo {

IR3 boris_push(
    const IR3& velocity, const double& tildeOref, const double& B, const IR3& b,
    const double& dt) {
  double T = std::tan(0.5 * tildeOref * dt * B), S = 2 * T / (1 + T * T);
  IR3 v_prime = velocity + T * cross_product(velocity, b);
  return velocity + S * cross_product(v_prime, b);
}

IR3 boris_push(
    const IR3& velocity, const double& tildeOref, const double& tildeEref,
    const IR3& E, const double& B, const IR3& b, const double& dt) {
  IR3 half_E_impulse = (0.5 * tildeEref * dt) * E;
  IR3 v_minus = velocity + half_E_impulse;
  IR3 v_plus = boris_push(v_minus, tildeOref, B, b, dt);
  return v_plus + half_E_impulse;
}

}  // end namespace gyronimo
