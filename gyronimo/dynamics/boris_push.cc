// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Manuel Assunção.

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
    const IR3& v_old, const double& Oref, const double& Bmag,
    const IR3& Bversor, const double& dt) {
  // step 1
  double T = std::tan(0.5 * Oref * dt * Bmag);
  double S = 2 * T / (1 + T * T);
  IR3 v_prime = v_old + T * cross_product(v_old, Bversor);
  IR3 v_new = v_old + S * cross_product(v_prime, Bversor);

  return v_new;
}

IR3 boris_push(
    const IR3& v_old, const double& Oref, const IR3& Efield, const double& Bmag,
    const IR3& Bversor, const double& dt) {
  // step 1
  IR3 half_E_impulse = (0.5 * Oref * dt) * Efield;
  IR3 v_minus = v_old + half_E_impulse;

  // step 2
  double T = std::tan(0.5 * Oref * dt * Bmag);
  double S = 2 * T / (1 + T * T);
  IR3 v_prime = v_minus + T * cross_product(v_minus, Bversor);
  IR3 v_plus = v_minus + S * cross_product(v_prime, Bversor);

  // step 3
  IR3 v_new = v_plus + half_E_impulse;

  return v_new;
}

}  // end namespace gyronimo