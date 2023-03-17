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

// @morphism_cartesian.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_cartesian.hh>

namespace gyronimo {

IR3 morphism_cartesian::operator()(const IR3& q) const { return q; }
IR3 morphism_cartesian::inverse(const IR3& x) const { return x; }
IR3 morphism_cartesian::translation(const IR3& q, const IR3& delta) const {
  return q + delta;
}
dIR3 morphism_cartesian::del(const IR3& q) const {
  return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
}
ddIR3 morphism_cartesian::ddel(const IR3& q) const {
  return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}
double morphism_cartesian::jacobian(const IR3& q) const { return 1.0; }
dIR3 morphism_cartesian::del_inverse(const IR3& q) const {
  return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
}
IR3 morphism_cartesian::to_covariant(const IR3& A, const IR3& q) const {
  return A;
}
IR3 morphism_cartesian::to_contravariant(const IR3& A, const IR3& q) const {
  return A;
}
IR3 morphism_cartesian::from_covariant(const IR3& A, const IR3& q) const {
  return A;
}
IR3 morphism_cartesian::from_contravariant(const IR3& A, const IR3& q) const {
  return A;
}

}  // end namespace gyronimo
