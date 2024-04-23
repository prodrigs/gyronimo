// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2024 Paulo Rodrigues and Manuel Assunção.

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
///@file

#include <gyronimo/core/contraction.hh>

namespace gyronimo {

//! Inner product @f$\mathbf{A}\cdot\mathbf{B} = A_i B^i = A^i B_i@f$.
double inner_product(const IR3& A, const IR3& B) {
  return A[IR3::u] * B[IR3::u] + A[IR3::v] * B[IR3::v] + A[IR3::w] * B[IR3::w];
}

//! Cartesian cross product.
/*!
    Implements the rule @f$(\mathbf{A} \times \mathbf{B})^k = (\mathbf{A} \times
    \mathbf{B})_k = \epsilon^{kij} A_i B_j = \epsilon_{kij} A^i B^j @f$, with
    @f$\epsilon@f$ the Levi-Civita symbol and valid for cartesian components
    only.
*/
IR3 cross_product(const IR3& A, const IR3& B) {
  return {
      A[IR3::v] * B[IR3::w] - A[IR3::w] * B[IR3::v],
      A[IR3::w] * B[IR3::u] - A[IR3::u] * B[IR3::w],
      A[IR3::u] * B[IR3::v] - A[IR3::v] * B[IR3::u]};
}

//! Covariant cross product, A and B interpreted as contravariant.
/*!
    Implements the rule @f$ (\mathbf{A} \times \mathbf{B})_k = \sqrt{g} \,
    \epsilon_{kij} A^i B^j @f$, with @f$\epsilon_{kij}@f$ the Levi-Civita
    symbol.
*/
template<>
IR3 cross_product<covariant>(const IR3& A, const IR3& B, double jacobian) {
  return cross_product(A, B) * jacobian;
}

//! Contravariant cross product, A and B interpreted as covariant.
/*!
    Implements the rule @f$(\mathbf{A} \times \mathbf{B})^k = \sqrt{g}^{-1}
    \epsilon^{kij} A_i B_j @f$, with @f$\epsilon^{kij}@f$ the Levi-Civita
    symbol.
*/
template<>
IR3 cross_product<contravariant>(const IR3& A, const IR3& B, double jacobian) {
  return cross_product(A, B) / jacobian;
}

//! Contraction of a `SM3` and a `IR3` vector.
/*!
    Returns the objects @f$ C_i = g_{ij} B^j @f$ or @f$ C^i = g^{ij} B_j @f$.
    Be sure to respect the variances of `g` and `B`.
*/
IR3 contraction(const SM3& g, const IR3& B) {
  return {
      g[SM3::uu] * B[IR3::u] + g[SM3::uv] * B[IR3::v] + g[SM3::uw] * B[IR3::w],
      g[SM3::uv] * B[IR3::u] + g[SM3::vv] * B[IR3::v] + g[SM3::vw] * B[IR3::w],
      g[SM3::uw] * B[IR3::u] + g[SM3::vw] * B[IR3::v] + g[SM3::ww] * B[IR3::w]};
}

//! First-index contraction of a `dIR3` object and a `IR3` vector.
/*!
    Returns the object @f$ C_i = B_j \partial_i A^j = B_j {A^j}_{,i} @f$ or @f$
    C_i = B^j \partial_i A_j = B^j {A_j}_{,i} @f$,
    ```
    C[IR3::i] = A[dIR3::ui]*B[IR3::u]
          + A[dIR3::vi]*B[IR3::v] + A[dIR3::wi]*B[IR3::w]
    ```
    with `i = u, v, w`. Be sure to respect the variances of `A` and `B`.
*/
template<> IR3 contraction<first>(const dIR3& A, const IR3& B) {
  return {
      A[dIR3::uu] * B[IR3::u] + A[dIR3::vu] * B[IR3::v] +
          A[dIR3::wu] * B[IR3::w],
      A[dIR3::uv] * B[IR3::u] + A[dIR3::vv] * B[IR3::v] +
          A[dIR3::wv] * B[IR3::w],
      A[dIR3::uw] * B[IR3::u] + A[dIR3::vw] * B[IR3::v] +
          A[dIR3::ww] * B[IR3::w]};
}

//! Second-index contraction of a `dIR3` object and a `IR3` vector.
/*!
    Returns the object @f$ C^i = B^j \partial_j A^i = B^j {A^i}_{,j} @f$ or @f$
    C_i = B^j \partial_j A_i = B^j {A_i}_{,j} @f$
    ```
    C[IR3::i] = A[dIR3::iu]*B[IR3::u]
      + A[dIR3::iv]*B[IR3::v] + A[dIR3::iw]*B[IR3::w]
    ```
    with `i = u, v, w`. Be sure to respect the variances of `A` and `B`.
*/
template<> IR3 contraction<second>(const dIR3& A, const IR3& B) {
  return {
      A[dIR3::uu] * B[IR3::u] + A[dIR3::uv] * B[IR3::v] +
          A[dIR3::uw] * B[IR3::w],
      A[dIR3::vu] * B[IR3::u] + A[dIR3::vv] * B[IR3::v] +
          A[dIR3::vw] * B[IR3::w],
      A[dIR3::wu] * B[IR3::u] + A[dIR3::wv] * B[IR3::v] +
          A[dIR3::ww] * B[IR3::w]};
}

//! First-index contraction of a `ddIR3` with a `SM3` matrix.
/*!
    Returns the object @f$ C^k_{,ij} = \partial_{ij} A_m B^{mk} = A_{m,ij}
    B^{mk} @f$ or @f$ C_{k,ij} = \partial_{ij} A^m B_{mk} = {A^m}_{,ij} B_{mk}
    @f$,
    ```
    C[ddIR3::ijk] = B[SM3::im] * A[ddIR3::mjk]
    ```
    where summation is implicit in the index `m`, with `m,i,j,k = u, v, w`.
    Replacements SM3::ij -> SM3::ji and ddIR3::ikj -> ddIR3::ijk have been
    explicitly performed. Be sure to respect the variances of `A` and `B`.
*/
template<> ddIR3 contraction<first>(const ddIR3& A, const SM3& B) {
  return {
      B[SM3::uu] * A[ddIR3::uuu] + B[SM3::uv] * A[ddIR3::vuu] +
          B[SM3::uw] * A[ddIR3::wuu],
      B[SM3::uu] * A[ddIR3::uuv] + B[SM3::uv] * A[ddIR3::vuv] +
          B[SM3::uw] * A[ddIR3::wuv],
      B[SM3::uu] * A[ddIR3::uuw] + B[SM3::uv] * A[ddIR3::vuw] +
          B[SM3::uw] * A[ddIR3::wuw],
      B[SM3::uu] * A[ddIR3::uvv] + B[SM3::uv] * A[ddIR3::vvv] +
          B[SM3::uw] * A[ddIR3::wvv],
      B[SM3::uu] * A[ddIR3::uvw] + B[SM3::uv] * A[ddIR3::vvw] +
          B[SM3::uw] * A[ddIR3::wvw],
      B[SM3::uu] * A[ddIR3::uww] + B[SM3::uv] * A[ddIR3::vww] +
          B[SM3::uw] * A[ddIR3::www],
      B[SM3::uv] * A[ddIR3::uuu] + B[SM3::vv] * A[ddIR3::vuu] +
          B[SM3::vw] * A[ddIR3::wuu],
      B[SM3::uv] * A[ddIR3::uuv] + B[SM3::vv] * A[ddIR3::vuv] +
          B[SM3::vw] * A[ddIR3::wuv],
      B[SM3::uv] * A[ddIR3::uuw] + B[SM3::vv] * A[ddIR3::vuw] +
          B[SM3::vw] * A[ddIR3::wuw],
      B[SM3::uv] * A[ddIR3::uvv] + B[SM3::vv] * A[ddIR3::vvv] +
          B[SM3::vw] * A[ddIR3::wvv],
      B[SM3::uv] * A[ddIR3::uvw] + B[SM3::vv] * A[ddIR3::vvw] +
          B[SM3::vw] * A[ddIR3::wvw],
      B[SM3::uv] * A[ddIR3::uww] + B[SM3::vv] * A[ddIR3::vww] +
          B[SM3::vw] * A[ddIR3::www],
      B[SM3::uw] * A[ddIR3::uuu] + B[SM3::vw] * A[ddIR3::vuu] +
          B[SM3::ww] * A[ddIR3::wuu],
      B[SM3::uw] * A[ddIR3::uuv] + B[SM3::vw] * A[ddIR3::vuv] +
          B[SM3::ww] * A[ddIR3::wuv],
      B[SM3::uw] * A[ddIR3::uuw] + B[SM3::vw] * A[ddIR3::vuw] +
          B[SM3::ww] * A[ddIR3::wuw],
      B[SM3::uw] * A[ddIR3::uvv] + B[SM3::vw] * A[ddIR3::vvv] +
          B[SM3::ww] * A[ddIR3::wvv],
      B[SM3::uw] * A[ddIR3::uvw] + B[SM3::vw] * A[ddIR3::vvw] +
          B[SM3::ww] * A[ddIR3::wvw],
      B[SM3::uw] * A[ddIR3::uww] + B[SM3::vw] * A[ddIR3::vww] +
          B[SM3::ww] * A[ddIR3::www]};
}

//! First-index contraction of a `dIR3` with a `ddIR3` objects.
/*!
    Returns the object @f$ C_{ijk} = {A^m}_{,i} B_{m,jk} @f$ or @f$ C_{kij} =
    {A_m}_{,i} {B^m}_{,jk} @f$,
    ```
    C[ddIR3::ijk] = A[dIR3::mi] * B[ddIR3::mjk]
    ```
    Be sure to respect the variances of `A` and `B`.
*/
template<> ddIR3 contraction<first>(const dIR3& A, const ddIR3& B) {
  return {
      A[dIR3::uu] * B[ddIR3::uuu] + A[dIR3::vu] * B[ddIR3::vuu] +
          A[dIR3::wu] * B[ddIR3::wuu],
      A[dIR3::uu] * B[ddIR3::uuv] + A[dIR3::vu] * B[ddIR3::vuv] +
          A[dIR3::wu] * B[ddIR3::wuv],
      A[dIR3::uu] * B[ddIR3::uuw] + A[dIR3::vu] * B[ddIR3::vuw] +
          A[dIR3::wu] * B[ddIR3::wuw],
      A[dIR3::uu] * B[ddIR3::uvv] + A[dIR3::vu] * B[ddIR3::vvv] +
          A[dIR3::wu] * B[ddIR3::wvv],
      A[dIR3::uu] * B[ddIR3::uvw] + A[dIR3::vu] * B[ddIR3::vvw] +
          A[dIR3::wu] * B[ddIR3::wvw],
      A[dIR3::uu] * B[ddIR3::uww] + A[dIR3::vu] * B[ddIR3::vww] +
          A[dIR3::wu] * B[ddIR3::www],
      A[dIR3::uv] * B[ddIR3::uuu] + A[dIR3::vv] * B[ddIR3::vuu] +
          A[dIR3::wv] * B[ddIR3::wuu],
      A[dIR3::uv] * B[ddIR3::uuv] + A[dIR3::vv] * B[ddIR3::vuv] +
          A[dIR3::wv] * B[ddIR3::wuv],
      A[dIR3::uv] * B[ddIR3::uuw] + A[dIR3::vv] * B[ddIR3::vuw] +
          A[dIR3::wv] * B[ddIR3::wuw],
      A[dIR3::uv] * B[ddIR3::uvv] + A[dIR3::vv] * B[ddIR3::vvv] +
          A[dIR3::wv] * B[ddIR3::wvv],
      A[dIR3::uv] * B[ddIR3::uvw] + A[dIR3::vv] * B[ddIR3::vvw] +
          A[dIR3::wv] * B[ddIR3::wvw],
      A[dIR3::uv] * B[ddIR3::uww] + A[dIR3::vv] * B[ddIR3::vww] +
          A[dIR3::wv] * B[ddIR3::www],
      A[dIR3::uw] * B[ddIR3::uuu] + A[dIR3::vw] * B[ddIR3::vuu] +
          A[dIR3::ww] * B[ddIR3::wuu],
      A[dIR3::uw] * B[ddIR3::uuv] + A[dIR3::vw] * B[ddIR3::vuv] +
          A[dIR3::ww] * B[ddIR3::wuv],
      A[dIR3::uw] * B[ddIR3::uuw] + A[dIR3::vw] * B[ddIR3::vuw] +
          A[dIR3::ww] * B[ddIR3::wuw],
      A[dIR3::uw] * B[ddIR3::uvv] + A[dIR3::vw] * B[ddIR3::vvv] +
          A[dIR3::ww] * B[ddIR3::wvv],
      A[dIR3::uw] * B[ddIR3::uvw] + A[dIR3::vw] * B[ddIR3::vvw] +
          A[dIR3::ww] * B[ddIR3::wvw],
      A[dIR3::uw] * B[ddIR3::uww] + A[dIR3::vw] * B[ddIR3::vww] +
          A[dIR3::ww] * B[ddIR3::www]};
}

//! Second-index contraction of a `dIR3` with a `ddIR3` objects.
/*!
    Returns the object @f$ {C^i}_{jk} = {A^i}_{,m} {B^m}_{,jk} @f$,
    ```
    C[ddIR3::ijk] = A[dIR3::im] * B[ddIR3::mjk]
    ```
    where summation is implicit in the index `m`, with `m,i,j,k = u, v, w`. Be
    sure to respect the variances of `A` and `B`.
*/
template<> ddIR3 contraction<second>(const dIR3& A, const ddIR3& B) {
  return {
      A[dIR3::uu] * B[ddIR3::uuu] + A[dIR3::uv] * B[ddIR3::vuu] +
          A[dIR3::uw] * B[ddIR3::wuu],
      A[dIR3::uu] * B[ddIR3::uuv] + A[dIR3::uv] * B[ddIR3::vuv] +
          A[dIR3::uw] * B[ddIR3::wuv],
      A[dIR3::uu] * B[ddIR3::uuw] + A[dIR3::uv] * B[ddIR3::vuw] +
          A[dIR3::uw] * B[ddIR3::wuw],
      A[dIR3::uu] * B[ddIR3::uvv] + A[dIR3::uv] * B[ddIR3::vvv] +
          A[dIR3::uw] * B[ddIR3::wvv],
      A[dIR3::uu] * B[ddIR3::uvw] + A[dIR3::uv] * B[ddIR3::vvw] +
          A[dIR3::uw] * B[ddIR3::wvw],
      A[dIR3::uu] * B[ddIR3::uww] + A[dIR3::uv] * B[ddIR3::vww] +
          A[dIR3::uw] * B[ddIR3::www],
      A[dIR3::vu] * B[ddIR3::uuu] + A[dIR3::vv] * B[ddIR3::vuu] +
          A[dIR3::vw] * B[ddIR3::wuu],
      A[dIR3::vu] * B[ddIR3::uuv] + A[dIR3::vv] * B[ddIR3::vuv] +
          A[dIR3::vw] * B[ddIR3::wuv],
      A[dIR3::vu] * B[ddIR3::uuw] + A[dIR3::vv] * B[ddIR3::vuw] +
          A[dIR3::vw] * B[ddIR3::wuw],
      A[dIR3::vu] * B[ddIR3::uvv] + A[dIR3::vv] * B[ddIR3::vvv] +
          A[dIR3::vw] * B[ddIR3::wvv],
      A[dIR3::vu] * B[ddIR3::uvw] + A[dIR3::vv] * B[ddIR3::vvw] +
          A[dIR3::vw] * B[ddIR3::wvw],
      A[dIR3::vu] * B[ddIR3::uww] + A[dIR3::vv] * B[ddIR3::vww] +
          A[dIR3::vw] * B[ddIR3::www],
      A[dIR3::wu] * B[ddIR3::uuu] + A[dIR3::wv] * B[ddIR3::vuu] +
          A[dIR3::ww] * B[ddIR3::wuu],
      A[dIR3::wu] * B[ddIR3::uuv] + A[dIR3::wv] * B[ddIR3::vuv] +
          A[dIR3::ww] * B[ddIR3::wuv],
      A[dIR3::wu] * B[ddIR3::uuw] + A[dIR3::wv] * B[ddIR3::vuw] +
          A[dIR3::ww] * B[ddIR3::wuw],
      A[dIR3::wu] * B[ddIR3::uvv] + A[dIR3::wv] * B[ddIR3::vvv] +
          A[dIR3::ww] * B[ddIR3::wvv],
      A[dIR3::wu] * B[ddIR3::uvw] + A[dIR3::wv] * B[ddIR3::vvw] +
          A[dIR3::ww] * B[ddIR3::wvw],
      A[dIR3::wu] * B[ddIR3::uww] + A[dIR3::wv] * B[ddIR3::vww] +
          A[dIR3::ww] * B[ddIR3::www]};
}

//! First-index contraction of a `dSM3` object and a `IR3` vector.
/*!
    Returns the object @f$ C_{ij} = \partial_j A_{ki} B^k = A_{ki,j} B^k @f$ or
    @f$ {C^i}_j = \partial_j A^{ki} B_k = {A^{ki}}_{,j} B_k @f$,
    ```
    C[dIR3::ij] = A[dSM3::uij]*B[IR3::u]
           + A[dSM3::vij]*B[IR3::v] + A[dSM3::wij]*B[IR3::w]
    ```
    with `i,j = u, v, w`. Be sure to respect the variances of `A` and `B`.
*/
template<> dIR3 contraction<first>(const dSM3& A, const IR3& B) {
  return {
      A[dSM3::uuu] * B[IR3::u] + A[dSM3::uvu] * B[IR3::v] +
          A[dSM3::uwu] * B[IR3::w],
      A[dSM3::uuv] * B[IR3::u] + A[dSM3::uvv] * B[IR3::v] +
          A[dSM3::uwv] * B[IR3::w],
      A[dSM3::uuw] * B[IR3::u] + A[dSM3::uvw] * B[IR3::v] +
          A[dSM3::uww] * B[IR3::w],
      A[dSM3::uvu] * B[IR3::u] + A[dSM3::vvu] * B[IR3::v] +
          A[dSM3::vwu] * B[IR3::w],
      A[dSM3::uvv] * B[IR3::u] + A[dSM3::vvv] * B[IR3::v] +
          A[dSM3::vwv] * B[IR3::w],
      A[dSM3::uvw] * B[IR3::u] + A[dSM3::vvw] * B[IR3::v] +
          A[dSM3::vww] * B[IR3::w],
      A[dSM3::uwu] * B[IR3::u] + A[dSM3::vwu] * B[IR3::v] +
          A[dSM3::wwu] * B[IR3::w],
      A[dSM3::uwv] * B[IR3::u] + A[dSM3::vwv] * B[IR3::v] +
          A[dSM3::wwv] * B[IR3::w],
      A[dSM3::uww] * B[IR3::u] + A[dSM3::vww] * B[IR3::v] +
          A[dSM3::www] * B[IR3::w]};
}

//! Second-index contraction of a `dSM3` object and a `IR3` vector.
template<> dIR3 contraction<second>(const dSM3& A, const IR3& B) {
  return contraction<first>(A, B);  // takes advantage of the dSM3 symmetry.
}

//! Third-index contraction of a `dSM3` object and a `IR3` vector.
/*!
    Returns the object @f$ C_{ij} = \partial_k A_{ij} B^k = A_{ij,k} B^k @f$ or
    @f$ C^{ij} = \partial_k A^{ij} B^k = {A^{ij}}_{,k} B_k @f$,
    ```
    C[dIR3::ij] = A[dSM3::iju]*B[IR3::u]
        + A[dSM3::ijv]*B[IR3::v] + A[dSM3::ijw]*B[IR3::w]
    ```
    with `i,j = u, v, w`. Be sure to respect the variances of `A` and `B`.
*/
template<> dIR3 contraction<third>(const dSM3& A, const IR3& B) {
  return {
      A[dSM3::uuu] * B[IR3::u] + A[dSM3::uuv] * B[IR3::v] +
          A[dSM3::uuw] * B[IR3::w],
      A[dSM3::uvu] * B[IR3::u] + A[dSM3::uvv] * B[IR3::v] +
          A[dSM3::uvw] * B[IR3::w],
      A[dSM3::uwu] * B[IR3::u] + A[dSM3::uwv] * B[IR3::v] +
          A[dSM3::uww] * B[IR3::w],
      A[dSM3::uvu] * B[IR3::u] + A[dSM3::uvv] * B[IR3::v] +
          A[dSM3::uvw] * B[IR3::w],
      A[dSM3::vvu] * B[IR3::u] + A[dSM3::vvv] * B[IR3::v] +
          A[dSM3::vvw] * B[IR3::w],
      A[dSM3::vwu] * B[IR3::u] + A[dSM3::vwv] * B[IR3::v] +
          A[dSM3::vww] * B[IR3::w],
      A[dSM3::uwu] * B[IR3::u] + A[dSM3::uwv] * B[IR3::v] +
          A[dSM3::uww] * B[IR3::w],
      A[dSM3::vwu] * B[IR3::u] + A[dSM3::vwv] * B[IR3::v] +
          A[dSM3::vww] * B[IR3::w],
      A[dSM3::wwu] * B[IR3::u] + A[dSM3::wwv] * B[IR3::v] +
          A[dSM3::www] * B[IR3::w]};
}

//! First-index contraction of a `dIR3` with a `dSM3` objects.
/*!
    Returns the object @f$ C_{ij} = \partial_j A^k B_{ki} = {A^k}_{,j} B_{ki}
    @f$ or @f$ {C^i}_j = \partial_j A_k B^{ki} = A_{k,j} B^{k,i} @f$,
    ```
    C[dIR3::ij] = B[SM3::iu]*A[dIR3::uj]
         + B[SM3::iv]*A[dIR3::vj] + B[SM3::iw]*A[dIR3::wj]
    ```
    with `i,j = u, v, w`. Be sure to respect the variances of `A` and `B`.
*/
template<> dIR3 contraction<first>(const dIR3& A, const SM3& B) {
  return {
      B[SM3::uu] * A[dIR3::uu] + B[SM3::uv] * A[dIR3::vu] +
          B[SM3::uw] * A[dIR3::wu],
      B[SM3::uu] * A[dIR3::uv] + B[SM3::uv] * A[dIR3::vv] +
          B[SM3::uw] * A[dIR3::wv],
      B[SM3::uu] * A[dIR3::uw] + B[SM3::uv] * A[dIR3::vw] +
          B[SM3::uw] * A[dIR3::ww],
      B[SM3::uv] * A[dIR3::uu] + B[SM3::vv] * A[dIR3::vu] +
          B[SM3::vw] * A[dIR3::wu],
      B[SM3::uv] * A[dIR3::uv] + B[SM3::vv] * A[dIR3::vv] +
          B[SM3::vw] * A[dIR3::wv],
      B[SM3::uv] * A[dIR3::uw] + B[SM3::vv] * A[dIR3::vw] +
          B[SM3::vw] * A[dIR3::ww],
      B[SM3::uw] * A[dIR3::uu] + B[SM3::vw] * A[dIR3::vu] +
          B[SM3::ww] * A[dIR3::wu],
      B[SM3::uw] * A[dIR3::uv] + B[SM3::vw] * A[dIR3::vv] +
          B[SM3::ww] * A[dIR3::wv],
      B[SM3::uw] * A[dIR3::uw] + B[SM3::vw] * A[dIR3::vw] +
          B[SM3::ww] * A[dIR3::ww]};
}

//! Second-index contraction of a `dIR3` with a `dSM3` objects.
/*!
    Returns the object @f$ C^{ij} = \partial_k A^j B^{ki} = {A^j}_{,k} B^{ki}
    @f$,
    ```
    C[dIR3::ij] = B[SM3::iu]*A[dIR3::ju]
         + B[SM3::iv]*A[dIR3::jv] + B[SM3::iw]*A[dIR3::jw]
    ```
    with `i,j = u, v, w`. Be sure to respect the variances of `A` and `B`.
*/
template<> dIR3 contraction<second>(const dIR3& A, const SM3& B) {
  return {
      B[SM3::uu] * A[dIR3::uu] + B[SM3::uv] * A[dIR3::uv] +
          B[SM3::uw] * A[dIR3::uw],
      B[SM3::uu] * A[dIR3::vu] + B[SM3::uv] * A[dIR3::vv] +
          B[SM3::uw] * A[dIR3::vw],
      B[SM3::uu] * A[dIR3::wu] + B[SM3::uv] * A[dIR3::wv] +
          B[SM3::uw] * A[dIR3::ww],
      B[SM3::uv] * A[dIR3::uu] + B[SM3::vv] * A[dIR3::uv] +
          B[SM3::vw] * A[dIR3::uw],
      B[SM3::uv] * A[dIR3::vu] + B[SM3::vv] * A[dIR3::vv] +
          B[SM3::vw] * A[dIR3::vw],
      B[SM3::uv] * A[dIR3::wu] + B[SM3::vv] * A[dIR3::wv] +
          B[SM3::vw] * A[dIR3::ww],
      B[SM3::uw] * A[dIR3::uu] + B[SM3::vw] * A[dIR3::uv] +
          B[SM3::ww] * A[dIR3::uw],
      B[SM3::uw] * A[dIR3::vu] + B[SM3::vw] * A[dIR3::vv] +
          B[SM3::ww] * A[dIR3::vw],
      B[SM3::uw] * A[dIR3::wu] + B[SM3::vw] * A[dIR3::wv] +
          B[SM3::ww] * A[dIR3::ww]};
}

}  // end namespace gyronimo.
