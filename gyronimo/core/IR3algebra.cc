// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @IR3algebra.cc, this file is part of ::gyronimo::

#include <gyronimo/core/IR3algebra.hh>

namespace gyronimo {

BinOpTree<IR3, IR3, std::plus<double>> const operator+(
    const IR3& x1, const IR3& x2) {
  return BinOpTree<IR3, IR3, std::plus<double>>(x1, x2);
}
BinOpTree<IR3, double, std::plus<double>> const operator+(
    const IR3& x1, const double& x2) {
  return BinOpTree<IR3, double, std::plus<double>>(x1, x2);
}
BinOpTree<double, IR3, std::plus<double>> const operator+(
    const double& x1, const IR3& x2) {
  return BinOpTree<double, IR3, std::plus<double>>(x1, x2);
}
BinOpTree<IR3, IR3, std::minus<double>> const operator-(
    const IR3& x1, const IR3& x2) {
  return BinOpTree<IR3, IR3, std::minus<double>>(x1, x2);
}
BinOpTree<IR3, double, std::minus<double>> const operator-(
    const IR3& x1, const double& x2) {
  return BinOpTree<IR3, double, std::minus<double>>(x1, x2);
}
BinOpTree<double, IR3, std::minus<double>> const operator-(
    const double& x1, const IR3& x2) {
  return BinOpTree<double, IR3, std::minus<double>>(x1, x2);
}
BinOpTree<IR3, IR3, std::multiplies<double>> const operator*(
    const IR3& x1, const IR3& x2) {
  return BinOpTree<IR3, IR3, std::multiplies<double>>(x1, x2);
}
BinOpTree<IR3, double, std::multiplies<double>> const operator*(
    const IR3& x1, const double& x2) {
  return BinOpTree<IR3, double, std::multiplies<double>>(x1, x2);
}
BinOpTree<double, IR3, std::multiplies<double>> const operator*(
    const double& x1, const IR3& x2) {
  return BinOpTree<double, IR3, std::multiplies<double>>(x1, x2);
}
BinOpTree<IR3, IR3, std::divides<double>> const operator/(
    const IR3& x1, const IR3& x2) {
  return BinOpTree<IR3, IR3, std::divides<double>>(x1, x2);
}
BinOpTree<IR3, double, std::divides<double>> const operator/(
    const IR3& x1, const double& x2) {
  return BinOpTree<IR3, double, std::divides<double>>(x1, x2);
}
BinOpTree<double, IR3, std::divides<double>> const operator/(
    const double& x1, const IR3& x2) {
  return BinOpTree<double, IR3, std::divides<double>>(x1, x2);
}

} // end namespace gyronimo.
