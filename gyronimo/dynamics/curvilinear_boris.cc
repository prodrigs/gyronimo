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

// @curvilinear_boris.cc, this file is part of ::gyronimo::

#include <gyronimo/dynamics/curvilinear_boris.hh>

namespace gyronimo {

curvilinear_boris::state curvilinear_boris::do_step(
    const state& s, const double& time, const double& dt) const {
  IR3 q = this->get_position(s);
  IR3 updated_v = this->cartesian_velocity_update(s, time, dt);
  IR3 dot_q_star = this->my_morphism()->to_contravariant(updated_v, q);
  IR3 q_half_step = q + (0.5 * this->Lref() * dt) * dot_q_star;
  IR3 dot_q_half_step =
      this->my_morphism()->to_contravariant(updated_v, q_half_step);
  IR3 updated_q = q + (this->Lref() * dt) * dot_q_half_step;
  return this->generate_state(updated_q, updated_v);
}

}  // end namespace gyronimo
