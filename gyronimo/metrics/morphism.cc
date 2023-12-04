// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Paulo Rodrigues and Manuel Assunção.

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

// @morphism.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Jacobian @f$ \mathbf{e}_u \cdot (\mathbf{e}_v \times \mathbf{e}_w) @f$.
double morphism::jacobian(const IR3& q) const {
  dIR3 e = del(q);
  return e[dIR3::uu] * (e[dIR3::vv] * e[dIR3::ww] - e[dIR3::vw] * e[dIR3::wv]) +
      e[dIR3::uv] * (e[dIR3::vw] * e[dIR3::wu] - e[dIR3::vu] * e[dIR3::ww]) +
      e[dIR3::uw] * (e[dIR3::vu] * e[dIR3::wv] - e[dIR3::vv] * e[dIR3::wu]);
}

}  // end namespace gyronimo
