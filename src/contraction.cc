// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @contraction.cc

#include <gyronimo/core/contraction.hh>

namespace gyronimo {

//! Inner product: one of A or B must be covariant and the other contravariant.
double inner_product(const IR3& A, const IR3& B) {
  return A[IR3::u]*B[IR3::u] + A[IR3::v]*B[IR3::v] + A[IR3::w]*B[IR3::w];
}

//! Contravariant cross product: A and B must be covariant.
IR3 cross_product(const IR3& A, const IR3& B, double jacobian) {
  return {
    (A[IR3::v]*B[IR3::w] - A[IR3::w]*B[IR3::v]) / jacobian,
    (A[IR3::w]*B[IR3::u] - A[IR3::u]*B[IR3::w]) / jacobian,
    (A[IR3::u]*B[IR3::v] - A[IR3::v]*B[IR3::u]) / jacobian};
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

//! Second-index contration of a `dIR3` object and a `IR3` vector.
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

//! First-index contration of a `dSM3` object and a `IR3` vector.
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

//! Second-index contration of a `dSM3` object and a `IR3` vector.
/*!
    Returns the object
    ```
    C[dIR3::ij] = dA[dSM3::iuj]*B[IR3::u]
        + dA[dSM3::ivj]*B[IR3::v] + dA[dSM3::iwj]*B[IR3::w]

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

//! Third-index contration of a `dSM3` object and a `IR3` vector.
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

//! Second-by-first-index contration of `dSM3` and `dIR3` objects.
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

//! Contration of a `dSM3` and two `SM3` objects.
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
