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

dIR3 inverse(const dIR3& m) {
  double ideterminant = 1.0/(
      m[dIR3::uu]*(m[dIR3::vv]*m[dIR3::ww] - m[dIR3::vw]*m[dIR3::wv]) -
      m[dIR3::uv]*(m[dIR3::vu]*m[dIR3::ww] - m[dIR3::vw]*m[dIR3::wu]) +
      m[dIR3::uw]*(m[dIR3::vu]*m[dIR3::wv] - m[dIR3::vv]*m[dIR3::wu]));
  return {
    ideterminant*(m[dIR3::vv]*m[dIR3::ww] - m[dIR3::vw]*m[dIR3::wv]),
    ideterminant*(m[dIR3::uw]*m[dIR3::wv] - m[dIR3::uv]*m[dIR3::ww]),
    ideterminant*(m[dIR3::uv]*m[dIR3::vw] - m[dIR3::uw]*m[dIR3::vv]),
    ideterminant*(m[dIR3::vw]*m[dIR3::wu] - m[dIR3::vu]*m[dIR3::ww]),
    ideterminant*(m[dIR3::uu]*m[dIR3::ww] - m[dIR3::uw]*m[dIR3::wu]),
    ideterminant*(m[dIR3::uw]*m[dIR3::vu] - m[dIR3::uu]*m[dIR3::vw]),
    ideterminant*(m[dIR3::vu]*m[dIR3::wv] - m[dIR3::vv]*m[dIR3::wu]),
    ideterminant*(m[dIR3::uv]*m[dIR3::wu] - m[dIR3::uu]*m[dIR3::wv]),
    ideterminant*(m[dIR3::uu]*m[dIR3::vv] - m[dIR3::uv]*m[dIR3::vu])};
}

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
