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

// @morphism_cylindrical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_cylindrical.hh>

#include <cmath>

namespace gyronimo {

IR3 morphism_cylindrical::operator()(const IR3& q) const {
  return {
      Lref_ * q[IR3::u] * std::cos(q[IR3::v]),
      Lref_ * q[IR3::u] * std::sin(q[IR3::v]), Lref_ * q[IR3::w]};
}
IR3 morphism_cylindrical::inverse(const IR3& x) const {
  return {
      iLref_ * std::sqrt(x[IR3::u] * x[IR3::u] + x[IR3::v] * x[IR3::v]),
      std::atan2(x[IR3::v], x[IR3::u]), iLref_ * x[IR3::w]};
}
dIR3 morphism_cylindrical::del(const IR3& q) const {
  double Lref_sn = Lref_ * std::sin(q[IR3::v]);
  double Lref_cn = Lref_ * std::cos(q[IR3::v]);
  double r = q[IR3::u];
  return {Lref_cn, -r * Lref_sn, 0, Lref_sn, r * Lref_cn, 0, 0, 0, Lref_};
}
ddIR3 morphism_cylindrical::ddel(const IR3& q) const {
  double r = q[IR3::u];
  double Lref_sn = Lref_ * std::sin(q[IR3::v]);
  double Lref_cn = Lref_ * std::cos(q[IR3::v]);
  return {0, -Lref_sn, 0, -r * Lref_cn, 0, 0,
          0, Lref_cn, 0, -r * Lref_sn, 0, 0, 0, 0, 0, 0, 0, 0};
}
double morphism_cylindrical::jacobian(const IR3& q) const {
  return Lref3_ * q[IR3::u];
}
dIR3 morphism_cylindrical::del_inverse(const IR3& q) const {
  double iLref_sn = iLref_ * std::sin(q[IR3::v]);
  double iLref_cn = iLref_ * std::cos(q[IR3::v]);
  double ir = 1 / q[IR3::u];
  return {
      iLref_cn, iLref_sn, 0, -iLref_sn * ir, iLref_cn * ir, 0, 0, 0, iLref_};
}

}  // end namespace gyronimo
