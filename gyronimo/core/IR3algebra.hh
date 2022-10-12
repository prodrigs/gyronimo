// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2022 Paulo Rodrigues.

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

// @IR3algebra.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_IR3ALGEBRA
#define GYRONIMO_IR3ALGEBRA

#include <array>
#include <functional>

namespace gyronimo {

//! Vector in @f$\mathbb{R}^3@f$ (i.e., IR3).
/*!
    *Non-aggregate* type, with implemented support for list initialization.
    Vectorised algebraic operations are also supported by overloading the
    corresponding operators (+, -, *, /, +=, -=, *=, /=) and temporary storage
    is avoided with expression templates. The index operator `B[IR3::i]` returns
    the value @f$B^{i}@f$ (for contravariant vectors) or @f$B_{i}@f$ (for
    covariant vectors) with `i = u, v, w`.

    @todo check if BinOpTree can be moved into a separate namespace; check if
    the initialisator can be removed from the constructor; check if the
    explicit loop and operator[] can be removed/abstracted.
*/
class IR3 : public  std::array<double, 3> {
 public:
  enum index {u = 0, v = 1, w = 2};

  IR3(const std::initializer_list<double>& list)
      : std::array<double, 3>::array() {
    std::copy(list.begin(), list.end(), this->begin());
  }
  template <typename ET>
  IR3(const ET& expr) {
    for(size_t i = 0;i < 3;i++) (*this)[i] = expr[i];
  }

  template <typename ET>
  IR3& operator+=(const ET& expr) {
    for(size_t i = 0;i < 3;i++) (*this)[i] += expr[i];
    return *this;
  }
  template <typename ET>
  IR3& operator-=(const ET& expr) {
    for(size_t i = 0;i < 3;i++) (*this)[i] -= expr[i];
    return *this;
  }
  template <typename ET>
  IR3& operator*=(const ET& expr) {
    for(size_t i = 0;i < 3;i++) (*this)[i] *= expr[i];
    return *this;
  }
  template <typename ET>
  IR3& operator/=(const ET& expr) {
    for(size_t i = 0;i < 3;i++) (*this)[i] /= expr[i];
    return *this;
  }
  IR3& operator+=(double x) {
    for(size_t i = 0;i < 3;i++) (*this)[i] += x;
    return *this;
  }
  IR3& operator-=(double x) {
    for(size_t i = 0;i < 3;i++) (*this)[i] -= x;
    return *this;
  }
  IR3& operator*=(double x) {
    for(size_t i = 0;i < 3;i++) (*this)[i] *= x;
    return *this;
  }
  IR3& operator/=(double x) {
    for(size_t i = 0;i < 3;i++) (*this)[i] /= x;
    return *this;
  }
};

//! Partial derivatives of a @f$\mathbb{R}^3@f$ vector.
/*!
    **Aggregate** type supporting list initialization. Unlike `IR3`, no
    vectorised algebraic operations are supported. The index operator
    `B[dIR3::ij]` returns the value @f$\partial_j B^i@f$ (for contravariant
    vectors) or @f$\partial_j B_i@f$ (for covariant vectors) with `i,j = u, v,
    w`.
*/
class dIR3 {
 public:
  enum index {
    uu = 0, uv = 1, uw = 2, vu = 3, vv = 4, vw = 5, wu = 6, wv = 7, ww = 8};

  double& operator[](index i) {return data_[i];}
  const double& operator[](index i) const {return data_[i];}
  template <typename T>
  dIR3& operator=(const T& expr) {
    for(size_t i = 0;i < 9;i++) data_[i] = expr[i];
    return *this;
  }

  std::array<double, 9> data_;
};

//! Double partial derivatives of a @f$\mathbb{R}^3@f$ vector.
/*!
    **Aggregate** type supporting list initialization. Unlike `IR3`, no
    vectorised algebraic operations are supported. The index operator
    `B[ddIR3::ijk]` returns the value @f$\partial_k\partial_j B^i@f$ (for contravariant
    vectors) or @f$\partial_k\partial_j B_i@f$ (for covariant vectors) with `i,j,k = u, v,
    w`.
*/
class ddIR3 {
 public:
  enum index {
    uuu = 0, uuv = 1, uuw = 2, uvv = 3, uvw = 4, uww = 5,
    vuu = 6, vuv = 7, vuw = 8, vvv = 9, vvw = 10, vww = 11,
    wuu = 12, wuv = 13, wuw = 14, wvv = 15, wvw = 16, www = 17,
	};

  double& operator[](index i) {return data_[i];}
  const double& operator[](index i) const {return data_[i];}
  template <typename T>
  ddIR3& operator=(const T& expr) {
    for(size_t i = 0;i < 18;i++) data_[i] = expr[i];
    return *this;
  }

  std::array<double, 18> data_;
};

//! Inverse of a dIR3 matrix.
dIR3 inverse(const dIR3& m);

//! Binary operation between two arbitrary types with `operator[i]`.
/*!
    Used to build at **compile time** Expression Templates representing an
    arbitrary combination of indexed objects via binary operations. The objects
    being used must have `operator[i]` implemented.
*/
template <typename T1, typename T2, class binop>
class BinOpTree {
 public:
  BinOpTree(const T1& x1, const T2& x2) : x1_(x1), x2_(x2) {}
  double operator[](size_t i) const {return binop()(x1_[i], x2_[i]);}

 private:
   const T1& x1_;
   const T2& x2_;
};

//! Binary operation between an arbitrary type with `operator[i]` and a double.
template <typename T, class binop>
class BinOpTree<T, double, binop> {
 public:
  BinOpTree(const T& x1, const double& x2) : x1_(x1), x2_(x2) {}
  double operator[](size_t i) const {return binop()(x1_[i], x2_);}

 private:
   const T& x1_;
   const double& x2_;
};

//! Binary operation between a double and an arbitrary type with `operator[i]`.
template <typename T, class binop>
class BinOpTree<double, T, binop> {
 public:
  BinOpTree(const double& x1, const T& x2) : x1_(x1), x2_(x2) {}
  double operator[](size_t i) const {return binop()(x1_, x2_[i]);}

 private:
   const double& x1_;
   const T& x2_;
};

// Overloads the algebraic operators within gyronimo:: in order to restrict the
// building of Expression Templates to IR3 objects only.
//
// 1.1 IR3 vectorised algebra, templated addition:
template <typename T1, typename T2, typename binop>
auto const operator+(const BinOpTree<T1, T2, binop>& x1, const IR3& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop>, IR3, std::plus<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator+(const IR3& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      IR3, BinOpTree<T1, T2, binop>, std::plus<double>>(x1, x2);
}
template <typename T1, typename T2,
    typename T3, typename T4, typename binop1, typename binop2>
auto const operator+(
    const BinOpTree<T1, T2, binop1>& x1, const BinOpTree<T3, T4, binop2>& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop1>,
      BinOpTree<T3, T4, binop2>, std::plus<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator+(const BinOpTree<T1, T2, binop>& x1, const double& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop>, double, std::plus<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator+(const double& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      double, BinOpTree<T1, T2, binop>, std::plus<double>>(x1, x2);
}

// 1.2 IR3 vectorised algebra, non-templated addition declaration:
// (operators do not support implicit instatiation).
BinOpTree<IR3, IR3, std::plus<double>> const
  operator+(const IR3& x1, const IR3& x2);
BinOpTree<IR3, double, std::plus<double>> const
  operator+(const IR3& x1, const double& x2);
BinOpTree<double, IR3, std::plus<double>> const
  operator+(const double& x1, const IR3& x2);

// 2.1 IR3 vectorised algebra, templated subtraction:
template <typename T1, typename T2, typename binop>
auto const operator-(const BinOpTree<T1, T2, binop>& x1, const IR3& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop>, IR3, std::minus<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator-(const IR3& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      IR3, BinOpTree<T1, T2, binop>, std::minus<double>>(x1, x2);
}
template <typename T1, typename T2,
    typename T3, typename T4, typename binop1, typename binop2>
auto const operator-(
    const BinOpTree<T1, T2, binop1>& x1, const BinOpTree<T3, T4, binop2>& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop1>,
      BinOpTree<T3, T4, binop2>, std::minus<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator-(const BinOpTree<T1, T2, binop>& x1, const double& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop>, double, std::minus<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator-(const double& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      double, BinOpTree<T1, T2, binop>, std::minus<double>>(x1, x2);
}

// 2.2 IR3 vectorised algebra, non-templated subtraction declaration:
// (operators do not support implicit instatiation).
BinOpTree<IR3, IR3, std::minus<double>> const
  operator-(const IR3& x1, const IR3& x2);
BinOpTree<IR3, double, std::minus<double>> const
  operator-(const IR3& x1, const double& x2);
BinOpTree<double, IR3, std::minus<double>> const
  operator-(const double& x1, const IR3& x2);

// 3.1 IR3 vectorised algebra, templated multiplication:
template <typename T1, typename T2, typename binop>
auto const operator*(const BinOpTree<T1, T2, binop>& x1, const IR3& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop>, IR3, std::multiplies<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator*(const IR3& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      IR3, BinOpTree<T1, T2, binop>, std::multiplies<double>>(x1, x2);
}
template <typename T1, typename T2,
    typename T3, typename T4, typename binop1, typename binop2>
auto const operator*(
    const BinOpTree<T1, T2, binop1>& x1, const BinOpTree<T3, T4, binop2>& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop1>,
      BinOpTree<T3, T4, binop2>, std::multiplies<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator*(const BinOpTree<T1, T2, binop>& x1, const double& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop>, double, std::multiplies<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator*(const double& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      double, BinOpTree<T1, T2, binop>, std::multiplies<double>>(x1, x2);
}

// 3.2 IR3 vectorised algebra, non-templated multiplication declaration:
// (operators do not support implicit instatiation).
BinOpTree<IR3, IR3, std::multiplies<double>> const
  operator*(const IR3& x1, const IR3& x2);
BinOpTree<IR3, double, std::multiplies<double>> const
  operator*(const IR3& x1, const double& x2);
BinOpTree<double, IR3, std::multiplies<double>> const
  operator*(const double& x1, const IR3& x2);

// 4.1 IR3 vectorised algebra, templated division:
template <typename T1, typename T2, typename binop>
auto const operator/(const BinOpTree<T1, T2, binop>& x1, const IR3& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop>, IR3, std::divides<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator/(const IR3& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      IR3, BinOpTree<T1, T2, binop>, std::divides<double>>(x1, x2);
}
template <typename T1, typename T2,
    typename T3, typename T4, typename binop1, typename binop2>
auto const operator/(
    const BinOpTree<T1, T2, binop1>& x1, const BinOpTree<T3, T4, binop2>& x2) {
  return BinOpTree<
      BinOpTree<T1, T2, binop1>,
      BinOpTree<T3, T4, binop2>, std::divides<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator/(const BinOpTree<T1, T2, binop>& x1, const double& x2) {
  return BinOpTree<
    BinOpTree<T1, T2, binop>, double, std::divides<double>>(x1, x2);
}
template <typename T1, typename T2, typename binop>
auto const operator/(const double& x1, const BinOpTree<T1, T2, binop>& x2) {
  return BinOpTree<
      double, BinOpTree<T1, T2, binop>, std::divides<double>>(x1, x2);
}

// 4.2 IR3 vectorised algebra, non-templated division declaration:
// (operators do not support implicit instatiation).
BinOpTree<IR3, IR3, std::divides<double>> const
  operator/(const IR3& x1, const IR3& x2);
BinOpTree<IR3, double, std::divides<double>> const
  operator/(const IR3& x1, const double& x2);
BinOpTree<double, IR3, std::divides<double>> const
  operator/(const double& x1, const IR3& x2);

} // end namespace gyronimo.

#endif // end GYRONIMO_IR3ALGEBRA
