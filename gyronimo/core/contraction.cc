// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2022 Paulo Rodrigues and Manuel Assunção.

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

// @contraction.cc, this file is part of ::gyronimo::

#include <gyronimo/core/contraction.hh>

namespace gyronimo {

//! Inner product: one of A or B must be covariant and the other contravariant.
double inner_product(const IR3& A, const IR3& B) {
  return A[IR3::u]*B[IR3::u] + A[IR3::v]*B[IR3::v] + A[IR3::w]*B[IR3::w];
}

//! Cartesian cross product: A and B must be cartesian.
IR3 cross_product(const IR3& A, const IR3& B) {
  return {
    A[IR3::v]*B[IR3::w] - A[IR3::w]*B[IR3::v],
    A[IR3::w]*B[IR3::u] - A[IR3::u]*B[IR3::w],
    A[IR3::u]*B[IR3::v] - A[IR3::v]*B[IR3::u]
  };
}

//! Covariant cross product: A and B must be contravariant.
template<>
IR3 cross_product<covariant>(const IR3& A, const IR3& B, double jacobian) {
  return {
    (A[IR3::v]*B[IR3::w] - A[IR3::w]*B[IR3::v]) * jacobian,
    (A[IR3::w]*B[IR3::u] - A[IR3::u]*B[IR3::w]) * jacobian,
    (A[IR3::u]*B[IR3::v] - A[IR3::v]*B[IR3::u]) * jacobian
  };
}

//! Contravariant cross product: A and B must be covariant.
template<>
IR3 cross_product<contravariant>(const IR3& A, const IR3& B, double jacobian) {
  double ijacobian = 1.0/jacobian;
  return {
    (A[IR3::v]*B[IR3::w] - A[IR3::w]*B[IR3::v]) * ijacobian,
    (A[IR3::w]*B[IR3::u] - A[IR3::u]*B[IR3::w]) * ijacobian,
    (A[IR3::u]*B[IR3::v] - A[IR3::v]*B[IR3::u]) * ijacobian
  };
}

//! Contraction of a symmetric 3x3 matrix and a IR^3 vector.
IR3 contraction(const SM3& g, const IR3& B) {
  return {
    g[SM3::uu]*B[IR3::u] + g[SM3::uv]*B[IR3::v] + g[SM3::uw]*B[IR3::w],
    g[SM3::uv]*B[IR3::u] + g[SM3::vv]*B[IR3::v] + g[SM3::vw]*B[IR3::w],
    g[SM3::uw]*B[IR3::u] + g[SM3::vw]*B[IR3::v] + g[SM3::ww]*B[IR3::w]};
}

//! First-index contraction of a `dIR3` object and a `IR3` vector.
/*!
    Returns the object
    ```
    C[IR3::i] = dA[dIR3::ui]*B[IR3::u]
          + dA[dIR3::vi]*B[IR3::v] + dA[dIR3::wi]*B[IR3::w]
    ```
    with `i = u, v, w`.
*/
template <>
IR3 contraction<first>(const dIR3& dA, const IR3& B) {
  return {
    dA[dIR3::uu]*B[IR3::u] + dA[dIR3::vu]*B[IR3::v] + dA[dIR3::wu]*B[IR3::w],
    dA[dIR3::uv]*B[IR3::u] + dA[dIR3::vv]*B[IR3::v] + dA[dIR3::wv]*B[IR3::w],
    dA[dIR3::uw]*B[IR3::u] + dA[dIR3::vw]*B[IR3::v] + dA[dIR3::ww]*B[IR3::w]};
}

//! Second-index contraction of a `dIR3` object and a `IR3` vector.
/*!
    Returns the object
    ```
    C[IR3::i] = dA[dIR3::iu]*B[IR3::u]
      + dA[dIR3::iv]*B[IR3::v] + dA[dIR3::iw]*B[IR3::w]
    ```
    with `i = u, v, w`.
*/
template <>
IR3 contraction<second>(const dIR3& dA, const IR3& B) {
  return {
    dA[dIR3::uu]*B[IR3::u] + dA[dIR3::uv]*B[IR3::v] + dA[dIR3::uw]*B[IR3::w],
    dA[dIR3::vu]*B[IR3::u] + dA[dIR3::vv]*B[IR3::v] + dA[dIR3::vw]*B[IR3::w],
    dA[dIR3::wu]*B[IR3::u] + dA[dIR3::wv]*B[IR3::v] + dA[dIR3::ww]*B[IR3::w]};
}

//! First and second index contraction of a `ddIR3` and two `IR3` vectors.
/*!
    Returns the object
    ```
    C[IR3::k] = ddA[ddIR3::ijk] * B[IR3::i] * C[IR3::j]
    ```
    where summation is implicit in the indices `i,j`, with `i,j,k = u, v, w`.
    Replacements ddIR3::ikj -> ddIR3::ijk have been explicitly performed.
*/
template<>
IR3 contraction<first, second>(const ddIR3& ddA, const IR3& B, const IR3& C) {
  return {
    ddA[ddIR3::uuu] * B[IR3::u] * C[IR3::u] +
    ddA[ddIR3::uuv] * B[IR3::u] * C[IR3::v] +
    ddA[ddIR3::uuw] * B[IR3::u] * C[IR3::w] +
    ddA[ddIR3::vuu] * B[IR3::v] * C[IR3::u] +
    ddA[ddIR3::vuv] * B[IR3::v] * C[IR3::v] +
    ddA[ddIR3::vuw] * B[IR3::v] * C[IR3::w] +
    ddA[ddIR3::wuu] * B[IR3::w] * C[IR3::u] +
    ddA[ddIR3::wuv] * B[IR3::w] * C[IR3::v] +
    ddA[ddIR3::wuw] * B[IR3::w] * C[IR3::w], // C[IR3::u]
    ddA[ddIR3::uuv] * B[IR3::u] * C[IR3::u] +
    ddA[ddIR3::uvv] * B[IR3::u] * C[IR3::v] +
    ddA[ddIR3::uvw] * B[IR3::u] * C[IR3::w] +
    ddA[ddIR3::vuv] * B[IR3::v] * C[IR3::u] +
    ddA[ddIR3::vvv] * B[IR3::v] * C[IR3::v] +
    ddA[ddIR3::vvw] * B[IR3::v] * C[IR3::w] +
    ddA[ddIR3::wuv] * B[IR3::w] * C[IR3::u] +
    ddA[ddIR3::wvv] * B[IR3::w] * C[IR3::v] +
    ddA[ddIR3::wvw] * B[IR3::w] * C[IR3::w], // C[IR3::v]
    ddA[ddIR3::uuw] * B[IR3::u] * C[IR3::u] +
    ddA[ddIR3::uvw] * B[IR3::u] * C[IR3::v] +
    ddA[ddIR3::uww] * B[IR3::u] * C[IR3::w] +
    ddA[ddIR3::vuw] * B[IR3::v] * C[IR3::u] +
    ddA[ddIR3::vvw] * B[IR3::v] * C[IR3::v] +
    ddA[ddIR3::vww] * B[IR3::v] * C[IR3::w] +
    ddA[ddIR3::wuw] * B[IR3::w] * C[IR3::u] +
    ddA[ddIR3::wvw] * B[IR3::w] * C[IR3::v] +
    ddA[ddIR3::www] * B[IR3::w] * C[IR3::w] // C[IR3::w]
  };
}

//! Second and third index contraction of a `ddIR3` and two `IR3` vectors.
/*!
    Returns the object
    ```
    C[IR3::i] = ddA[ddIR3::ijk] * B[IR3::j] * C[IR3::k]
    ```
    where summation is implicit in the indices `j,k`, with `i,j,k = u, v, w`.
    Replacements ddIR3::ikj -> ddIR3::ijk have been explicitly performed.
*/
template<>
IR3 contraction<second, third>(const ddIR3& ddA, const IR3& B, const IR3& C) {
  return {
    ddA[ddIR3::uuu] *  B[IR3::u] * C[IR3::u] +
    ddA[ddIR3::uuv] * (B[IR3::u] * C[IR3::v] + B[IR3::v] * C[IR3::u]) +
    ddA[ddIR3::uuw] * (B[IR3::u] * C[IR3::w] + B[IR3::w] * C[IR3::u]) +
    ddA[ddIR3::uvv] *  B[IR3::v] * C[IR3::v] +
    ddA[ddIR3::uvw] * (B[IR3::v] * C[IR3::w] + B[IR3::w] * C[IR3::v]) +
    ddA[ddIR3::uww] *  B[IR3::w] * C[IR3::w], // C[IR3::u]
    ddA[ddIR3::vuu] *  B[IR3::u] * C[IR3::u] +
    ddA[ddIR3::vuv] * (B[IR3::u] * C[IR3::v] + B[IR3::v] * C[IR3::u]) +
    ddA[ddIR3::vuw] * (B[IR3::u] * C[IR3::w] + B[IR3::w] * C[IR3::u]) +
    ddA[ddIR3::vvv] *  B[IR3::v] * C[IR3::v] +
    ddA[ddIR3::vvw] * (B[IR3::v] * C[IR3::w] + B[IR3::w] * C[IR3::v]) +
    ddA[ddIR3::vww] *  B[IR3::w] * C[IR3::w], // C[IR3::v]
    ddA[ddIR3::wuu] *  B[IR3::u] * C[IR3::u] +
    ddA[ddIR3::wuv] * (B[IR3::u] * C[IR3::v] + B[IR3::v] * C[IR3::u]) +
    ddA[ddIR3::wuw] * (B[IR3::u] * C[IR3::w] + B[IR3::w] * C[IR3::u]) +
    ddA[ddIR3::wvv] *  B[IR3::v] * C[IR3::v] +
    ddA[ddIR3::wvw] * (B[IR3::v] * C[IR3::w] + B[IR3::w] * C[IR3::v]) +
    ddA[ddIR3::www] *  B[IR3::w] * C[IR3::w] // C[IR3::w]
  };
}

//! First-index contraction of a `ddIR3` by second index of a `SM3` matrix.
/*!
    Returns the object
    ```
    ddC[ddIR3::ijk] = g[SM3::im] * ddA[ddIR3::mjk]
    ```
    where summation is implicit in the index `m`, with `m,i,j,k = u, v, w`.
    Replacements SM3::ij -> SM3::ji and ddIR3::ikj -> ddIR3::ijk have been
    explicitly performed.
*/
template<>
ddIR3 contraction<second, first>(const SM3& g, const ddIR3& ddA) {
  return {
    g[SM3::uu]*ddA[ddIR3::uuu] + g[SM3::uv]*ddA[ddIR3::vuu] + g[SM3::uw]*ddA[ddIR3::wuu], // uuu
    g[SM3::uu]*ddA[ddIR3::uuv] + g[SM3::uv]*ddA[ddIR3::vuv] + g[SM3::uw]*ddA[ddIR3::wuv], // uuv
    g[SM3::uu]*ddA[ddIR3::uuw] + g[SM3::uv]*ddA[ddIR3::vuw] + g[SM3::uw]*ddA[ddIR3::wuw], // uuw
    g[SM3::uu]*ddA[ddIR3::uvv] + g[SM3::uv]*ddA[ddIR3::vvv] + g[SM3::uw]*ddA[ddIR3::wvv], // uvv
    g[SM3::uu]*ddA[ddIR3::uvw] + g[SM3::uv]*ddA[ddIR3::vvw] + g[SM3::uw]*ddA[ddIR3::wvw], // uvw
    g[SM3::uu]*ddA[ddIR3::uww] + g[SM3::uv]*ddA[ddIR3::vww] + g[SM3::uw]*ddA[ddIR3::www], // uww
    g[SM3::uv]*ddA[ddIR3::uuu] + g[SM3::vv]*ddA[ddIR3::vuu] + g[SM3::vw]*ddA[ddIR3::wuu], // vuu
    g[SM3::uv]*ddA[ddIR3::uuv] + g[SM3::vv]*ddA[ddIR3::vuv] + g[SM3::vw]*ddA[ddIR3::wuv], // vuv
    g[SM3::uv]*ddA[ddIR3::uuw] + g[SM3::vv]*ddA[ddIR3::vuw] + g[SM3::vw]*ddA[ddIR3::wuw], // vuw
    g[SM3::uv]*ddA[ddIR3::uvv] + g[SM3::vv]*ddA[ddIR3::vvv] + g[SM3::vw]*ddA[ddIR3::wvv], // vvv
    g[SM3::uv]*ddA[ddIR3::uvw] + g[SM3::vv]*ddA[ddIR3::vvw] + g[SM3::vw]*ddA[ddIR3::wvw], // vvw
    g[SM3::uv]*ddA[ddIR3::uww] + g[SM3::vv]*ddA[ddIR3::vww] + g[SM3::vw]*ddA[ddIR3::www], // vww
    g[SM3::uw]*ddA[ddIR3::uuu] + g[SM3::vw]*ddA[ddIR3::vuu] + g[SM3::ww]*ddA[ddIR3::wuu], // wuu
    g[SM3::uw]*ddA[ddIR3::uuv] + g[SM3::vw]*ddA[ddIR3::vuv] + g[SM3::ww]*ddA[ddIR3::wuv], // wuv
    g[SM3::uw]*ddA[ddIR3::uuw] + g[SM3::vw]*ddA[ddIR3::vuw] + g[SM3::ww]*ddA[ddIR3::wuw], // wuw
    g[SM3::uw]*ddA[ddIR3::uvv] + g[SM3::vw]*ddA[ddIR3::vvv] + g[SM3::ww]*ddA[ddIR3::wvv], // wvv
    g[SM3::uw]*ddA[ddIR3::uvw] + g[SM3::vw]*ddA[ddIR3::vvw] + g[SM3::ww]*ddA[ddIR3::wvw], // wvw
    g[SM3::uw]*ddA[ddIR3::uww] + g[SM3::vw]*ddA[ddIR3::vww] + g[SM3::ww]*ddA[ddIR3::www], // www
  };
}

//! First-by-first index contraction of `dIR3` and `ddIR3` objects.
/*!
    Returns the object
    ```
    ddC[ddIR3::ijk] = dA[dIR3::mi] * ddB[ddIR3::mjk]
    ```
    where summation is implicit in the index `m`, with `m,i,j,k = u, v, w`.
    Replacements ddIR3::ikj -> ddIR3::ijk have been explicitly performed.
*/
template<>
ddIR3 contraction<first, first>(const dIR3& dA, const ddIR3& ddB) {
  return {
    dA[dIR3::uu]*ddB[ddIR3::uuu] + dA[dIR3::vu]*ddB[ddIR3::vuu] + dA[dIR3::wu]*ddB[ddIR3::wuu], // uuu
    dA[dIR3::uu]*ddB[ddIR3::uuv] + dA[dIR3::vu]*ddB[ddIR3::vuv] + dA[dIR3::wu]*ddB[ddIR3::wuv], // uuv
    dA[dIR3::uu]*ddB[ddIR3::uuw] + dA[dIR3::vu]*ddB[ddIR3::vuw] + dA[dIR3::wu]*ddB[ddIR3::wuw], // uuw
    dA[dIR3::uu]*ddB[ddIR3::uvv] + dA[dIR3::vu]*ddB[ddIR3::vvv] + dA[dIR3::wu]*ddB[ddIR3::wvv], // uvv
    dA[dIR3::uu]*ddB[ddIR3::uvw] + dA[dIR3::vu]*ddB[ddIR3::vvw] + dA[dIR3::wu]*ddB[ddIR3::wvw], // uvw
    dA[dIR3::uu]*ddB[ddIR3::uww] + dA[dIR3::vu]*ddB[ddIR3::vww] + dA[dIR3::wu]*ddB[ddIR3::www], // uww
    dA[dIR3::uv]*ddB[ddIR3::uuu] + dA[dIR3::vv]*ddB[ddIR3::vuu] + dA[dIR3::wv]*ddB[ddIR3::wuu], // vuu
    dA[dIR3::uv]*ddB[ddIR3::uuv] + dA[dIR3::vv]*ddB[ddIR3::vuv] + dA[dIR3::wv]*ddB[ddIR3::wuv], // vuv
    dA[dIR3::uv]*ddB[ddIR3::uuw] + dA[dIR3::vv]*ddB[ddIR3::vuw] + dA[dIR3::wv]*ddB[ddIR3::wuw], // vuw
    dA[dIR3::uv]*ddB[ddIR3::uvv] + dA[dIR3::vv]*ddB[ddIR3::vvv] + dA[dIR3::wv]*ddB[ddIR3::wvv], // vvv
    dA[dIR3::uv]*ddB[ddIR3::uvw] + dA[dIR3::vv]*ddB[ddIR3::vvw] + dA[dIR3::wv]*ddB[ddIR3::wvw], // vvw
    dA[dIR3::uv]*ddB[ddIR3::uww] + dA[dIR3::vv]*ddB[ddIR3::vww] + dA[dIR3::wv]*ddB[ddIR3::www], // vww
    dA[dIR3::uw]*ddB[ddIR3::uuu] + dA[dIR3::vw]*ddB[ddIR3::vuu] + dA[dIR3::ww]*ddB[ddIR3::wuu], // wuu
    dA[dIR3::uw]*ddB[ddIR3::uuv] + dA[dIR3::vw]*ddB[ddIR3::vuv] + dA[dIR3::ww]*ddB[ddIR3::wuv], // wuv
    dA[dIR3::uw]*ddB[ddIR3::uuw] + dA[dIR3::vw]*ddB[ddIR3::vuw] + dA[dIR3::ww]*ddB[ddIR3::wuw], // wuw
    dA[dIR3::uw]*ddB[ddIR3::uvv] + dA[dIR3::vw]*ddB[ddIR3::vvv] + dA[dIR3::ww]*ddB[ddIR3::wvv], // wvv
    dA[dIR3::uw]*ddB[ddIR3::uvw] + dA[dIR3::vw]*ddB[ddIR3::vvw] + dA[dIR3::ww]*ddB[ddIR3::wvw], // wvw
    dA[dIR3::uw]*ddB[ddIR3::uww] + dA[dIR3::vw]*ddB[ddIR3::vww] + dA[dIR3::ww]*ddB[ddIR3::www], // www
  };
}

//! Second-by-first index contraction of `dIR3` and `ddIR3` objects.
/*!
    Returns the object
    ```
    ddC[ddIR3::ijk] = dA[dIR3::im] * ddB[ddIR3::mjk]
    ```
    where summation is implicit in the index `m`, with `m,i,j,k = u, v, w`.
    Replacements ddIR3::ikj -> ddIR3::ijk have been explicitly performed.
*/
template<>
ddIR3 contraction<second, first>(const dIR3& dA, const ddIR3& ddB) {
  return {
    dA[dIR3::uu]*ddB[ddIR3::uuu] + dA[dIR3::uv]*ddB[ddIR3::vuu] + dA[dIR3::uw]*ddB[ddIR3::wuu], // uuu
    dA[dIR3::uu]*ddB[ddIR3::uuv] + dA[dIR3::uv]*ddB[ddIR3::vuv] + dA[dIR3::uw]*ddB[ddIR3::wuv], // uuv
    dA[dIR3::uu]*ddB[ddIR3::uuw] + dA[dIR3::uv]*ddB[ddIR3::vuw] + dA[dIR3::uw]*ddB[ddIR3::wuw], // uuw
    dA[dIR3::uu]*ddB[ddIR3::uvv] + dA[dIR3::uv]*ddB[ddIR3::vvv] + dA[dIR3::uw]*ddB[ddIR3::wvv], // uvv
    dA[dIR3::uu]*ddB[ddIR3::uvw] + dA[dIR3::uv]*ddB[ddIR3::vvw] + dA[dIR3::uw]*ddB[ddIR3::wvw], // uvw
    dA[dIR3::uu]*ddB[ddIR3::uww] + dA[dIR3::uv]*ddB[ddIR3::vww] + dA[dIR3::uw]*ddB[ddIR3::www], // uww
    dA[dIR3::vu]*ddB[ddIR3::uuu] + dA[dIR3::vv]*ddB[ddIR3::vuu] + dA[dIR3::vw]*ddB[ddIR3::wuu], // vuu
    dA[dIR3::vu]*ddB[ddIR3::uuv] + dA[dIR3::vv]*ddB[ddIR3::vuv] + dA[dIR3::vw]*ddB[ddIR3::wuv], // vuv
    dA[dIR3::vu]*ddB[ddIR3::uuw] + dA[dIR3::vv]*ddB[ddIR3::vuw] + dA[dIR3::vw]*ddB[ddIR3::wuw], // vuw
    dA[dIR3::vu]*ddB[ddIR3::uvv] + dA[dIR3::vv]*ddB[ddIR3::vvv] + dA[dIR3::vw]*ddB[ddIR3::wvv], // vvv
    dA[dIR3::vu]*ddB[ddIR3::uvw] + dA[dIR3::vv]*ddB[ddIR3::vvw] + dA[dIR3::vw]*ddB[ddIR3::wvw], // vvw
    dA[dIR3::vu]*ddB[ddIR3::uww] + dA[dIR3::vv]*ddB[ddIR3::vww] + dA[dIR3::vw]*ddB[ddIR3::www], // vww
    dA[dIR3::wu]*ddB[ddIR3::uuu] + dA[dIR3::wv]*ddB[ddIR3::vuu] + dA[dIR3::ww]*ddB[ddIR3::wuu], // wuu
    dA[dIR3::wu]*ddB[ddIR3::uuv] + dA[dIR3::wv]*ddB[ddIR3::vuv] + dA[dIR3::ww]*ddB[ddIR3::wuv], // wuv
    dA[dIR3::wu]*ddB[ddIR3::uuw] + dA[dIR3::wv]*ddB[ddIR3::vuw] + dA[dIR3::ww]*ddB[ddIR3::wuw], // wuw
    dA[dIR3::wu]*ddB[ddIR3::uvv] + dA[dIR3::wv]*ddB[ddIR3::vvv] + dA[dIR3::ww]*ddB[ddIR3::wvv], // wvv
    dA[dIR3::wu]*ddB[ddIR3::uvw] + dA[dIR3::wv]*ddB[ddIR3::vvw] + dA[dIR3::ww]*ddB[ddIR3::wvw], // wvw
    dA[dIR3::wu]*ddB[ddIR3::uww] + dA[dIR3::wv]*ddB[ddIR3::vww] + dA[dIR3::ww]*ddB[ddIR3::www], // www
  };
}

//! First-index contraction of a `dSM3` object and a `IR3` vector.
/*!
    Returns the object
    ```
    C[dIR3::ij] = dA[dSM3::uij]*B[IR3::u]
           + dA[dSM3::vij]*B[IR3::v] + dA[dSM3::wij]*B[IR3::w]
    ```
    with `i,j = u, v, w`. Replacements dSM3::ijk -> dSM3::jik have been
    explicitly performed.
*/
template<>
dIR3 contraction<first>(const dSM3& dA, const IR3& B) {
  return {
    dA[dSM3::uuu]*B[IR3::u] + dA[dSM3::uvu]*B[IR3::v] + dA[dSM3::uwu]*B[IR3::w],
    dA[dSM3::uuv]*B[IR3::u] + dA[dSM3::uvv]*B[IR3::v] + dA[dSM3::uwv]*B[IR3::w],
    dA[dSM3::uuw]*B[IR3::u] + dA[dSM3::uvw]*B[IR3::v] + dA[dSM3::uww]*B[IR3::w],
    dA[dSM3::uvu]*B[IR3::u] + dA[dSM3::vvu]*B[IR3::v] + dA[dSM3::vwu]*B[IR3::w],
    dA[dSM3::uvv]*B[IR3::u] + dA[dSM3::vvv]*B[IR3::v] + dA[dSM3::vwv]*B[IR3::w],
    dA[dSM3::uvw]*B[IR3::u] + dA[dSM3::vvw]*B[IR3::v] + dA[dSM3::vww]*B[IR3::w],
    dA[dSM3::uwu]*B[IR3::u] + dA[dSM3::vwu]*B[IR3::v] + dA[dSM3::wwu]*B[IR3::w],
    dA[dSM3::uwv]*B[IR3::u] + dA[dSM3::vwv]*B[IR3::v] + dA[dSM3::wwv]*B[IR3::w],
    dA[dSM3::uww]*B[IR3::u] + dA[dSM3::vww]*B[IR3::v] + dA[dSM3::www]*B[IR3::w]
  };
}

//! Second-index contraction of a `dSM3` object and a `IR3` vector.
/*!
    Returns the object
    ```
    C[dIR3::ij] = dA[dSM3::iuj]*B[IR3::u]
        + dA[dSM3::ivj]*B[IR3::v] + dA[dSM3::iwj]*B[IR3::w]
    ```
    with `i,j = u, v, w`. Replacements dSM3::ijk -> dSM3::jik have been
    explicitly performed.

    @todo check if this <second> can be replaced by <first> because dSM3 are
    symmetric in their first and second indices.
*/
template<>
dIR3 contraction<second>(const dSM3& dA, const IR3& B) {
  return {
    dA[dSM3::uuu]*B[IR3::u] + dA[dSM3::uvu]*B[IR3::v] + dA[dSM3::uwu]*B[IR3::w],
    dA[dSM3::uuv]*B[IR3::u] + dA[dSM3::uvv]*B[IR3::v] + dA[dSM3::uwv]*B[IR3::w],
    dA[dSM3::uuw]*B[IR3::u] + dA[dSM3::uvw]*B[IR3::v] + dA[dSM3::uww]*B[IR3::w],
    dA[dSM3::uvu]*B[IR3::u] + dA[dSM3::vvu]*B[IR3::v] + dA[dSM3::vwu]*B[IR3::w],
    dA[dSM3::uvv]*B[IR3::u] + dA[dSM3::vvv]*B[IR3::v] + dA[dSM3::vwv]*B[IR3::w],
    dA[dSM3::uvw]*B[IR3::u] + dA[dSM3::vvw]*B[IR3::v] + dA[dSM3::vww]*B[IR3::w],
    dA[dSM3::uwu]*B[IR3::u] + dA[dSM3::vwu]*B[IR3::v] + dA[dSM3::wwu]*B[IR3::w],
    dA[dSM3::uwv]*B[IR3::u] + dA[dSM3::vwv]*B[IR3::v] + dA[dSM3::wwv]*B[IR3::w],
    dA[dSM3::uww]*B[IR3::u] + dA[dSM3::vww]*B[IR3::v] + dA[dSM3::www]*B[IR3::w]
  };
}

//! Third-index contraction of a `dSM3` object and a `IR3` vector.
/*!
    Returns the object
    ```
    C[dIR3::ij] = dA[dSM3::iju]*B[IR3::u]
        + dA[dSM3::ijv]*B[IR3::v] + dA[dSM3::ijw]*B[IR3::w]
    ```
    with `i,j = u, v, w`. Replacements dSM3::ijk -> dSM3::jik have been
    explicitly performed.
*/
template<>
dIR3 contraction<third>(const dSM3& dA, const IR3& B) {
  return {
    dA[dSM3::uuu]*B[IR3::u] + dA[dSM3::uuv]*B[IR3::v] + dA[dSM3::uuw]*B[IR3::w],
    dA[dSM3::uvu]*B[IR3::u] + dA[dSM3::uvv]*B[IR3::v] + dA[dSM3::uvw]*B[IR3::w],
    dA[dSM3::uwu]*B[IR3::u] + dA[dSM3::uwv]*B[IR3::v] + dA[dSM3::uww]*B[IR3::w],
    dA[dSM3::uvu]*B[IR3::u] + dA[dSM3::uvv]*B[IR3::v] + dA[dSM3::uvw]*B[IR3::w],
    dA[dSM3::vvu]*B[IR3::u] + dA[dSM3::vvv]*B[IR3::v] + dA[dSM3::vvw]*B[IR3::w],
    dA[dSM3::vwu]*B[IR3::u] + dA[dSM3::vwv]*B[IR3::v] + dA[dSM3::vww]*B[IR3::w],
    dA[dSM3::uwu]*B[IR3::u] + dA[dSM3::uwv]*B[IR3::v] + dA[dSM3::uww]*B[IR3::w],
    dA[dSM3::vwu]*B[IR3::u] + dA[dSM3::vwv]*B[IR3::v] + dA[dSM3::vww]*B[IR3::w],
    dA[dSM3::wwu]*B[IR3::u] + dA[dSM3::wwv]*B[IR3::v] + dA[dSM3::www]*B[IR3::w]
  };
}

//! Second-by-first-index contraction of `dSM3` and `dIR3` objects.
/*!
    Returns the object
    ```
    C[dIR3::ij] = g[SM3::iu]*dB[dIR3::uj]
         + g[SM3::iv]*dB[dIR3::vj] + g[SM3::iw]*dB[dIR3::wj]
    ```
    with `i,j = u, v, w`. Replacements dSM3::ijk -> dSM3::jik have been
    explicitly performed.
*/
template<>
dIR3 contraction<second, first>(const SM3& g, const dIR3& dB) {
  return {
    g[SM3::uu]*dB[dIR3::uu] + g[SM3::uv]*dB[dIR3::vu] + g[SM3::uw]*dB[dIR3::wu],
    g[SM3::uu]*dB[dIR3::uv] + g[SM3::uv]*dB[dIR3::vv] + g[SM3::uw]*dB[dIR3::wv],
    g[SM3::uu]*dB[dIR3::uw] + g[SM3::uv]*dB[dIR3::vw] + g[SM3::uw]*dB[dIR3::ww],
    g[SM3::uv]*dB[dIR3::uu] + g[SM3::vv]*dB[dIR3::vu] + g[SM3::vw]*dB[dIR3::wu],
    g[SM3::uv]*dB[dIR3::uv] + g[SM3::vv]*dB[dIR3::vv] + g[SM3::vw]*dB[dIR3::wv],
    g[SM3::uv]*dB[dIR3::uw] + g[SM3::vv]*dB[dIR3::vw] + g[SM3::vw]*dB[dIR3::ww],
    g[SM3::uw]*dB[dIR3::uu] + g[SM3::vw]*dB[dIR3::vu] + g[SM3::ww]*dB[dIR3::wu],
    g[SM3::uw]*dB[dIR3::uv] + g[SM3::vw]*dB[dIR3::vv] + g[SM3::ww]*dB[dIR3::wv],
    g[SM3::uw]*dB[dIR3::uw] + g[SM3::vw]*dB[dIR3::vw] + g[SM3::ww]*dB[dIR3::ww]
  };
}

//! Contraction of a `dSM3` and two `SM3` objects.
/*!
    Returns the object
    ```
    C[dSM3::ijk] = d[dSM3::mnk]*g[SM3::mi]*h[SM3::nj]
    ```
    where summation is implicit in the indices `m,n`, with `i,j = u, v, w`.
    Replacements dSM3::jik -> dSM3::ijk and SM3:ji -> SM3::ij have been
    explicitly performed.
*/
dSM3 contraction(const SM3& g, const dSM3& d, const SM3& h) {
  return {
    g[SM3::uu]*( // these are dSM3::uu?
      d[dSM3::uuu]*h[SM3::uu]+d[dSM3::uvu]*h[SM3::uv]+d[dSM3::uwu]*h[SM3::uw]) +
    g[SM3::uv]*(
      d[dSM3::uvu]*h[SM3::uu]+d[dSM3::vvu]*h[SM3::uv]+d[dSM3::vwu]*h[SM3::uw]) +
    g[SM3::uw]*(
      d[dSM3::uwu]*h[SM3::uu]+d[dSM3::vwu]*h[SM3::uv]+d[dSM3::wwu]*h[SM3::uw]),
    g[SM3::uu]*(
      d[dSM3::uuv]*h[SM3::uu]+d[dSM3::uvv]*h[SM3::uv]+d[dSM3::uwv]*h[SM3::uw]) +
    g[SM3::uv]*(
      d[dSM3::uvv]*h[SM3::uu]+d[dSM3::vvv]*h[SM3::uv]+d[dSM3::vwv]*h[SM3::uw]) +
    g[SM3::uw]*(
      d[dSM3::uwv]*h[SM3::uu]+d[dSM3::vwv]*h[SM3::uv]+d[dSM3::wwv]*h[SM3::uw]),
    g[SM3::uu]*(
      d[dSM3::uuw]*h[SM3::uu]+d[dSM3::uvw]*h[SM3::uv]+d[dSM3::uww]*h[SM3::uw]) +
    g[SM3::uv]*(
      d[dSM3::uvw]*h[SM3::uu]+d[dSM3::vvw]*h[SM3::uv]+d[dSM3::vww]*h[SM3::uw]) +
    g[SM3::uw]*(
      d[dSM3::uww]*h[SM3::uu]+d[dSM3::vww]*h[SM3::uv]+d[dSM3::www]*h[SM3::uw]),
    g[SM3::uu]*( // these are dSM3::uv?
      d[dSM3::uuu]*h[SM3::uv]+d[dSM3::uvu]*h[SM3::vv]+d[dSM3::uwu]*h[SM3::vw]) +
    g[SM3::uv]*(
      d[dSM3::uvu]*h[SM3::uv]+d[dSM3::vvu]*h[SM3::vv]+d[dSM3::vwu]*h[SM3::vw]) +
    g[SM3::uw]*(
      d[dSM3::uwu]*h[SM3::uv]+d[dSM3::vwu]*h[SM3::vv]+d[dSM3::wwu]*h[SM3::vw]),
    g[SM3::uu]*(
      d[dSM3::uuv]*h[SM3::uv]+d[dSM3::uvv]*h[SM3::vv]+d[dSM3::uwv]*h[SM3::vw]) +
    g[SM3::uv]*(
      d[dSM3::uvv]*h[SM3::uv]+d[dSM3::vvv]*h[SM3::vv]+d[dSM3::vwv]*h[SM3::vw]) +
    g[SM3::uw]*(
      d[dSM3::uwv]*h[SM3::uv]+d[dSM3::vwv]*h[SM3::vv]+d[dSM3::wwv]*h[SM3::vw]),
    g[SM3::uu]*(
      d[dSM3::uuw]*h[SM3::uv]+d[dSM3::uvw]*h[SM3::vv]+d[dSM3::uww]*h[SM3::vw]) +
    g[SM3::uv]*(
      d[dSM3::uvw]*h[SM3::uv]+d[dSM3::vvw]*h[SM3::vv]+d[dSM3::vww]*h[SM3::vw]) +
    g[SM3::uw]*(
      d[dSM3::uww]*h[SM3::uv]+d[dSM3::vww]*h[SM3::vv]+d[dSM3::www]*h[SM3::vw]),
    g[SM3::uu]*( // these are dSM3::uw?
      d[dSM3::uuu]*h[SM3::uw]+d[dSM3::uvu]*h[SM3::vw]+d[dSM3::uwu]*h[SM3::ww]) +
    g[SM3::uv]*(
      d[dSM3::uvu]*h[SM3::uw]+d[dSM3::vvu]*h[SM3::vw]+d[dSM3::vwu]*h[SM3::ww]) +
    g[SM3::uw]*(
      d[dSM3::uwu]*h[SM3::uw]+d[dSM3::vwu]*h[SM3::vw]+d[dSM3::wwu]*h[SM3::ww]),
    g[SM3::uu]*(
      d[dSM3::uuv]*h[SM3::uw]+d[dSM3::uvv]*h[SM3::vw]+d[dSM3::uwv]*h[SM3::ww]) +
    g[SM3::uv]*(
      d[dSM3::uvv]*h[SM3::uw]+d[dSM3::vvv]*h[SM3::vw]+d[dSM3::vwv]*h[SM3::ww]) +
    g[SM3::uw]*(
      d[dSM3::uwv]*h[SM3::uw]+d[dSM3::vwv]*h[SM3::vw]+d[dSM3::wwv]*h[SM3::ww]),
    g[SM3::uu]*(
      d[dSM3::uuw]*h[SM3::uw]+d[dSM3::uvw]*h[SM3::vw]+d[dSM3::uww]*h[SM3::ww]) +
    g[SM3::uv]*(
      d[dSM3::uvw]*h[SM3::uw]+d[dSM3::vvw]*h[SM3::vw]+d[dSM3::vww]*h[SM3::ww]) +
    g[SM3::uw]*(
      d[dSM3::uww]*h[SM3::uw]+d[dSM3::vww]*h[SM3::vw]+d[dSM3::www]*h[SM3::ww]),
    g[SM3::uv]*( // these are dSM3::vv?
      d[dSM3::uuu]*h[SM3::uv]+d[dSM3::uvu]*h[SM3::vv]+d[dSM3::uwu]*h[SM3::vw]) +
    g[SM3::vv]*(
      d[dSM3::uvu]*h[SM3::uv]+d[dSM3::vvu]*h[SM3::vv]+d[dSM3::vwu]*h[SM3::vw]) +
    g[SM3::vw]*(
      d[dSM3::uwu]*h[SM3::uv]+d[dSM3::vwu]*h[SM3::vv]+d[dSM3::wwu]*h[SM3::vw]),
    g[SM3::uv]*(
      d[dSM3::uuv]*h[SM3::uv]+d[dSM3::uvv]*h[SM3::vv]+d[dSM3::uwv]*h[SM3::vw]) +
    g[SM3::vv]*(
      d[dSM3::uvv]*h[SM3::uv]+d[dSM3::vvv]*h[SM3::vv]+d[dSM3::vwv]*h[SM3::vw]) +
    g[SM3::vw]*(
      d[dSM3::uwv]*h[SM3::uv]+d[dSM3::vwv]*h[SM3::vv]+d[dSM3::wwv]*h[SM3::vw]),
    g[SM3::uv]*(
      d[dSM3::uuw]*h[SM3::uv]+d[dSM3::uvw]*h[SM3::vv]+d[dSM3::uww]*h[SM3::vw]) +
    g[SM3::vv]*(
      d[dSM3::uvw]*h[SM3::uv]+d[dSM3::vvw]*h[SM3::vv]+d[dSM3::vww]*h[SM3::vw]) +
    g[SM3::vw]*(
      d[dSM3::uww]*h[SM3::uv]+d[dSM3::vww]*h[SM3::vv]+d[dSM3::www]*h[SM3::vw]),
    g[SM3::uv]*( // these are dSM3::vw?
      d[dSM3::uuu]*h[SM3::uw]+d[dSM3::uvu]*h[SM3::vw]+d[dSM3::uwu]*h[SM3::ww]) +
    g[SM3::vv]*(
      d[dSM3::uvu]*h[SM3::uw]+d[dSM3::vvu]*h[SM3::vw]+d[dSM3::vwu]*h[SM3::ww]) +
    g[SM3::vw]*(
      d[dSM3::uwu]*h[SM3::uw]+d[dSM3::vwu]*h[SM3::vw]+d[dSM3::wwu]*h[SM3::ww]),
    g[SM3::uv]*(
      d[dSM3::uuv]*h[SM3::uw]+d[dSM3::uvv]*h[SM3::vw]+d[dSM3::uwv]*h[SM3::ww]) +
    g[SM3::vv]*(
      d[dSM3::uvv]*h[SM3::uw]+d[dSM3::vvv]*h[SM3::vw]+d[dSM3::vwv]*h[SM3::ww]) +
    g[SM3::vw]*(
      d[dSM3::uwv]*h[SM3::uw]+d[dSM3::vwv]*h[SM3::vw]+d[dSM3::wwv]*h[SM3::ww]),
    g[SM3::uv]*(
      d[dSM3::uuw]*h[SM3::uw]+d[dSM3::uvw]*h[SM3::vw]+d[dSM3::uww]*h[SM3::ww]) +
    g[SM3::vv]*(
      d[dSM3::uvw]*h[SM3::uw]+d[dSM3::vvw]*h[SM3::vw]+d[dSM3::vww]*h[SM3::ww]) +
    g[SM3::vw]*(
      d[dSM3::uww]*h[SM3::uw]+d[dSM3::vww]*h[SM3::vw]+d[dSM3::www]*h[SM3::ww]),
    g[SM3::uw]*( // these are dSM3::ww?
      d[dSM3::uuu]*h[SM3::uw]+d[dSM3::uvu]*h[SM3::vw]+d[dSM3::uwu]*h[SM3::ww]) +
    g[SM3::vw]*(
      d[dSM3::uvu]*h[SM3::uw]+d[dSM3::vvu]*h[SM3::vw]+d[dSM3::vwu]*h[SM3::ww]) +
    g[SM3::ww]*(
      d[dSM3::uwu]*h[SM3::uw]+d[dSM3::vwu]*h[SM3::vw]+d[dSM3::wwu]*h[SM3::ww]),
    g[SM3::uw]*(
      d[dSM3::uuv]*h[SM3::uw]+d[dSM3::uvv]*h[SM3::vw]+d[dSM3::uwv]*h[SM3::ww]) +
    g[SM3::vw]*(
      d[dSM3::uvv]*h[SM3::uw]+d[dSM3::vvv]*h[SM3::vw]+d[dSM3::vwv]*h[SM3::ww]) +
    g[SM3::ww]*(
      d[dSM3::uwv]*h[SM3::uw]+d[dSM3::vwv]*h[SM3::vw]+d[dSM3::wwv]*h[SM3::ww]),
    g[SM3::uw]*(
      d[dSM3::uuw]*h[SM3::uw]+d[dSM3::uvw]*h[SM3::vw]+d[dSM3::uww]*h[SM3::ww]) +
    g[SM3::vw]*(
      d[dSM3::uvw]*h[SM3::uw]+d[dSM3::vvw]*h[SM3::vw]+d[dSM3::vww]*h[SM3::ww]) +
    g[SM3::ww]*(
      d[dSM3::uww]*h[SM3::uw]+d[dSM3::vww]*h[SM3::vw]+d[dSM3::www]*h[SM3::ww])};
}

} // end namespace gyronimo.
