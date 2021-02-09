// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @IR3algebra.cc

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
