// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues and Manuel Assunção.

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

// @metric_covariant.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_covariant.hh>

#include <cmath>

namespace gyronimo {

//! General-purpose implementation of the Jacobian.
double metric_covariant::jacobian(const IR3& r) const {
  SM3 g = (*this)(r);
  double det_g = g[SM3::uu]*g[SM3::vv]*g[SM3::ww] +
      2*g[SM3::uv]*g[SM3::uw]*g[SM3::vw] - g[SM3::uv]*g[SM3::uv]*g[SM3::ww] -
      g[SM3::uu]*g[SM3::vw]*g[SM3::vw] - g[SM3::uw]*g[SM3::uw]*g[SM3::vv];
  return std::sqrt(det_g);
}

//! General-purpose implementation of the Jacobian gradient.
IR3 metric_covariant::del_jacobian(const IR3& r) const {
  SM3 g = (*this)(r);
  dSM3 dg = this->del(r);
  double det_g = g[SM3::uu]*g[SM3::vv]*g[SM3::ww] +
      2*g[SM3::uv]*g[SM3::uw]*g[SM3::vw] - g[SM3::uv]*g[SM3::uv]*g[SM3::ww] -
      g[SM3::uu]*g[SM3::vw]*g[SM3::vw] - g[SM3::uw]*g[SM3::uw]*g[SM3::vv];
  IR3 grad_det_g = {
      2*g[SM3::vw]*(
          g[SM3::uw]*dg[dSM3::uvu] +
          g[SM3::uv]*dg[dSM3::uwu] - g[SM3::uu]*dg[dSM3::vwu]) -
      g[SM3::uv]*g[SM3::uv]*dg[dSM3::wwu] + g[SM3::vv]*(
          g[SM3::ww]*dg[dSM3::uuu] -
          2*g[SM3::uw]*dg[dSM3::uwu] + g[SM3::uu]*dg[dSM3::wwu]) -
      g[SM3::vw]*g[SM3::vw]*dg[dSM3::uuu] -
      2*g[SM3::uv]*g[SM3::ww]*dg[dSM3::uvu] - 
      g[SM3::uw]*g[SM3::uw]*dg[dSM3::vvu] +
      g[SM3::uu]*g[SM3::ww]*dg[dSM3::vvu] + 
      2*g[SM3::uv]*g[SM3::uw]*dg[dSM3::vwu],
      2*g[SM3::vw]*(
          g[SM3::uw]*dg[dSM3::uvv] +
          g[SM3::uv]*dg[dSM3::uwv] - g[SM3::uu]*dg[dSM3::vwv]) -
      g[SM3::uv]*g[SM3::uv]*dg[dSM3::wwv] + g[SM3::vv]*(
          g[SM3::ww]*dg[dSM3::uuv] -
          2*g[SM3::uw]*dg[dSM3::uwv] + g[SM3::uu]*dg[dSM3::wwv]) -
      g[SM3::vw]*g[SM3::vw]*dg[dSM3::uuv] -
      2*g[SM3::uv]*g[SM3::ww]*dg[dSM3::uvv] - 
      g[SM3::uw]*g[SM3::uw]*dg[dSM3::vvv] +
      g[SM3::uu]*g[SM3::ww]*dg[dSM3::vvv] + 
      2*g[SM3::uv]*g[SM3::uw]*dg[dSM3::vwv],
      2*g[SM3::vw]*(
          g[SM3::uw]*dg[dSM3::uvw] +
          g[SM3::uv]*dg[dSM3::uww] - g[SM3::uu]*dg[dSM3::vww]) -
      g[SM3::uv]*g[SM3::uv]*dg[dSM3::www] + g[SM3::vv]*(
          g[SM3::ww]*dg[dSM3::uuw] -
          2*g[SM3::uw]*dg[dSM3::uww] + g[SM3::uu]*dg[dSM3::www]) -
      g[SM3::vw]*g[SM3::vw]*dg[dSM3::uuw] -
      2*g[SM3::uv]*g[SM3::ww]*dg[dSM3::uvw] - 
      g[SM3::uw]*g[SM3::uw]*dg[dSM3::vvw] +
      g[SM3::uu]*g[SM3::ww]*dg[dSM3::vvw] + 
      2*g[SM3::uv]*g[SM3::uw]*dg[dSM3::vww]};
  return (0.5/std::sqrt(det_g))*grad_det_g;
}

//! General derivative @f$\partial_k g^{ij} = -g^{im}\partial_k g_{mn}g^{nj}@f$.
dSM3 metric_covariant::del_inverse(const IR3& q) const {
  dSM3 dg = this->del(q);
  SM3 g = this->inverse(q);
  SM3 g2 = {g[SM3::uu] * g[SM3::uu], g[SM3::uv] * g[SM3::uv],
      g[SM3::uw] * g[SM3::uw], g[SM3::vv] * g[SM3::vv], g[SM3::vw] * g[SM3::vw],
      g[SM3::ww] * g[SM3::ww]};
  return {
      -2 * dg[dSM3::uvu] * g[SM3::uu] * g[SM3::uv] -
          2 * dg[dSM3::uwu] * g[SM3::uu] * g[SM3::uw] -
          2 * dg[dSM3::vwu] * g[SM3::uv] * g[SM3::uw] -
          dg[dSM3::uuu] * g2[SM3::uu] - dg[dSM3::vvu] * g2[SM3::uv] -
          dg[dSM3::wwu] * g2[SM3::uw],
      -2 * dg[dSM3::uvv] * g[SM3::uu] * g[SM3::uv] -
          2 * dg[dSM3::uwv] * g[SM3::uu] * g[SM3::uw] -
          2 * dg[dSM3::vwv] * g[SM3::uv] * g[SM3::uw] -
          dg[dSM3::uuv] * g2[SM3::uu] - dg[dSM3::vvv] * g2[SM3::uv] -
          dg[dSM3::wwv] * g2[SM3::uw],
      -2 * dg[dSM3::uvw] * g[SM3::uu] * g[SM3::uv] -
          2 * dg[dSM3::uww] * g[SM3::uu] * g[SM3::uw] -
          2 * dg[dSM3::vww] * g[SM3::uv] * g[SM3::uw] -
          dg[dSM3::uuw] * g2[SM3::uu] - dg[dSM3::vvw] * g2[SM3::uv] -
          dg[dSM3::www] * g2[SM3::uw],
      -(g[SM3::uu] *
        (dg[dSM3::uuu] * g[SM3::uv] + dg[dSM3::uvu] * g[SM3::vv] +
         dg[dSM3::uwu] * g[SM3::vw])) -
          g[SM3::uv] *
              (dg[dSM3::uvu] * g[SM3::uv] + dg[dSM3::vvu] * g[SM3::vv] +
               dg[dSM3::vwu] * g[SM3::vw]) -
          g[SM3::uw] *
              (dg[dSM3::uwu] * g[SM3::uv] + dg[dSM3::vwu] * g[SM3::vv] +
               dg[dSM3::wwu] * g[SM3::vw]),
      -(g[SM3::uu] *
        (dg[dSM3::uuv] * g[SM3::uv] + dg[dSM3::uvv] * g[SM3::vv] +
         dg[dSM3::uwv] * g[SM3::vw])) -
          g[SM3::uv] *
              (dg[dSM3::uvv] * g[SM3::uv] + dg[dSM3::vvv] * g[SM3::vv] +
               dg[dSM3::vwv] * g[SM3::vw]) -
          g[SM3::uw] *
              (dg[dSM3::uwv] * g[SM3::uv] + dg[dSM3::vwv] * g[SM3::vv] +
               dg[dSM3::wwv] * g[SM3::vw]),
      -(g[SM3::uu] *
        (dg[dSM3::uuw] * g[SM3::uv] + dg[dSM3::uvw] * g[SM3::vv] +
         dg[dSM3::uww] * g[SM3::vw])) -
          g[SM3::uv] *
              (dg[dSM3::uvw] * g[SM3::uv] + dg[dSM3::vvw] * g[SM3::vv] +
               dg[dSM3::vww] * g[SM3::vw]) -
          g[SM3::uw] *
              (dg[dSM3::uww] * g[SM3::uv] + dg[dSM3::vww] * g[SM3::vv] +
               dg[dSM3::www] * g[SM3::vw]),
      -(g[SM3::uu] *
        (dg[dSM3::uuu] * g[SM3::uw] + dg[dSM3::uvu] * g[SM3::vw] +
         dg[dSM3::uwu] * g[SM3::ww])) -
          g[SM3::uv] *
              (dg[dSM3::uvu] * g[SM3::uw] + dg[dSM3::vvu] * g[SM3::vw] +
               dg[dSM3::vwu] * g[SM3::ww]) -
          g[SM3::uw] *
              (dg[dSM3::uwu] * g[SM3::uw] + dg[dSM3::vwu] * g[SM3::vw] +
               dg[dSM3::wwu] * g[SM3::ww]),
      -(g[SM3::uu] *
        (dg[dSM3::uuv] * g[SM3::uw] + dg[dSM3::uvv] * g[SM3::vw] +
         dg[dSM3::uwv] * g[SM3::ww])) -
          g[SM3::uv] *
              (dg[dSM3::uvv] * g[SM3::uw] + dg[dSM3::vvv] * g[SM3::vw] +
               dg[dSM3::vwv] * g[SM3::ww]) -
          g[SM3::uw] *
              (dg[dSM3::uwv] * g[SM3::uw] + dg[dSM3::vwv] * g[SM3::vw] +
               dg[dSM3::wwv] * g[SM3::ww]),
      -(g[SM3::uu] *
        (dg[dSM3::uuw] * g[SM3::uw] + dg[dSM3::uvw] * g[SM3::vw] +
         dg[dSM3::uww] * g[SM3::ww])) -
          g[SM3::uv] *
              (dg[dSM3::uvw] * g[SM3::uw] + dg[dSM3::vvw] * g[SM3::vw] +
               dg[dSM3::vww] * g[SM3::ww]) -
          g[SM3::uw] *
              (dg[dSM3::uww] * g[SM3::uw] + dg[dSM3::vww] * g[SM3::vw] +
               dg[dSM3::www] * g[SM3::ww]),
      -2 * dg[dSM3::uvu] * g[SM3::uv] * g[SM3::vv] -
          2 * dg[dSM3::uwu] * g[SM3::uv] * g[SM3::vw] -
          2 * dg[dSM3::vwu] * g[SM3::vv] * g[SM3::vw] -
          dg[dSM3::uuu] * g2[SM3::uv] - dg[dSM3::vvu] * g2[SM3::vv] -
          dg[dSM3::wwu] * g2[SM3::vw],
      -2 * dg[dSM3::uvv] * g[SM3::uv] * g[SM3::vv] -
          2 * dg[dSM3::uwv] * g[SM3::uv] * g[SM3::vw] -
          2 * dg[dSM3::vwv] * g[SM3::vv] * g[SM3::vw] -
          dg[dSM3::uuv] * g2[SM3::uv] - dg[dSM3::vvv] * g2[SM3::vv] -
          dg[dSM3::wwv] * g2[SM3::vw],
      -2 * dg[dSM3::uvw] * g[SM3::uv] * g[SM3::vv] -
          2 * dg[dSM3::uww] * g[SM3::uv] * g[SM3::vw] -
          2 * dg[dSM3::vww] * g[SM3::vv] * g[SM3::vw] -
          dg[dSM3::uuw] * g2[SM3::uv] - dg[dSM3::vvw] * g2[SM3::vv] -
          dg[dSM3::www] * g2[SM3::vw],
      -(g[SM3::uv] *
        (dg[dSM3::uuu] * g[SM3::uw] + dg[dSM3::uvu] * g[SM3::vw] +
         dg[dSM3::uwu] * g[SM3::ww])) -
          g[SM3::vv] *
              (dg[dSM3::uvu] * g[SM3::uw] + dg[dSM3::vvu] * g[SM3::vw] +
               dg[dSM3::vwu] * g[SM3::ww]) -
          g[SM3::vw] *
              (dg[dSM3::uwu] * g[SM3::uw] + dg[dSM3::vwu] * g[SM3::vw] +
               dg[dSM3::wwu] * g[SM3::ww]),
      -(g[SM3::uv] *
        (dg[dSM3::uuv] * g[SM3::uw] + dg[dSM3::uvv] * g[SM3::vw] +
         dg[dSM3::uwv] * g[SM3::ww])) -
          g[SM3::vv] *
              (dg[dSM3::uvv] * g[SM3::uw] + dg[dSM3::vvv] * g[SM3::vw] +
               dg[dSM3::vwv] * g[SM3::ww]) -
          g[SM3::vw] *
              (dg[dSM3::uwv] * g[SM3::uw] + dg[dSM3::vwv] * g[SM3::vw] +
               dg[dSM3::wwv] * g[SM3::ww]),
      -(g[SM3::uv] *
        (dg[dSM3::uuw] * g[SM3::uw] + dg[dSM3::uvw] * g[SM3::vw] +
         dg[dSM3::uww] * g[SM3::ww])) -
          g[SM3::vv] *
              (dg[dSM3::uvw] * g[SM3::uw] + dg[dSM3::vvw] * g[SM3::vw] +
               dg[dSM3::vww] * g[SM3::ww]) -
          g[SM3::vw] *
              (dg[dSM3::uww] * g[SM3::uw] + dg[dSM3::vww] * g[SM3::vw] +
               dg[dSM3::www] * g[SM3::ww]),
      -2 * dg[dSM3::uvu] * g[SM3::uw] * g[SM3::vw] -
          2 * dg[dSM3::uwu] * g[SM3::uw] * g[SM3::ww] -
          2 * dg[dSM3::vwu] * g[SM3::vw] * g[SM3::ww] -
          dg[dSM3::uuu] * g2[SM3::uw] - dg[dSM3::vvu] * g2[SM3::vw] -
          dg[dSM3::wwu] * g2[SM3::ww],
      -2 * dg[dSM3::uvv] * g[SM3::uw] * g[SM3::vw] -
          2 * dg[dSM3::uwv] * g[SM3::uw] * g[SM3::ww] -
          2 * dg[dSM3::vwv] * g[SM3::vw] * g[SM3::ww] -
          dg[dSM3::uuv] * g2[SM3::uw] - dg[dSM3::vvv] * g2[SM3::vw] -
          dg[dSM3::wwv] * g2[SM3::ww],
      -2 * dg[dSM3::uvw] * g[SM3::uw] * g[SM3::vw] -
          2 * dg[dSM3::uww] * g[SM3::uw] * g[SM3::ww] -
          2 * dg[dSM3::vww] * g[SM3::vw] * g[SM3::ww] -
          dg[dSM3::uuw] * g2[SM3::uw] - dg[dSM3::vvw] * g2[SM3::vw] -
          dg[dSM3::www] * g2[SM3::ww]};
}

//! General Christoffel symbol @f$\Gamma_{ijk}@f$.
/*!
    Implements the rule @f$ \Gamma_{ijk} = \frac{1}{2} ( \partial_k g_{ij} +
    \partial_j g_{ik} - \partial_i g_{jk}) @f$.
*/
ddIR3 metric_covariant::christoffel_first_kind(const IR3& q) const {
  dSM3 dg = del(q);
  return {
      0.5 * (dg[dSM3::uuu] + dg[dSM3::uuu] - dg[dSM3::uuu]),  // uuu
      0.5 * (dg[dSM3::uuv] + dg[dSM3::uvu] - dg[dSM3::uvu]),  // uuv
      0.5 * (dg[dSM3::uuw] + dg[dSM3::uwu] - dg[dSM3::uwu]),  // uuw
      0.5 * (dg[dSM3::uvv] + dg[dSM3::uvv] - dg[dSM3::vvu]),  // uvv
      0.5 * (dg[dSM3::uvw] + dg[dSM3::uwv] - dg[dSM3::vwu]),  // uvw
      0.5 * (dg[dSM3::uww] + dg[dSM3::uww] - dg[dSM3::wwu]),  // uww
      0.5 * (dg[dSM3::uvu] + dg[dSM3::uvu] - dg[dSM3::uuv]),  // vuu
      0.5 * (dg[dSM3::uvv] + dg[dSM3::vvu] - dg[dSM3::uvv]),  // vuv
      0.5 * (dg[dSM3::uvw] + dg[dSM3::vwu] - dg[dSM3::uwv]),  // vuw
      0.5 * (dg[dSM3::vvv] + dg[dSM3::vvv] - dg[dSM3::vvv]),  // vvv
      0.5 * (dg[dSM3::vvw] + dg[dSM3::vwv] - dg[dSM3::vwv]),  // vvw
      0.5 * (dg[dSM3::vww] + dg[dSM3::vww] - dg[dSM3::wwv]),  // vww
      0.5 * (dg[dSM3::uwu] + dg[dSM3::uwu] - dg[dSM3::uuw]),  // wuu
      0.5 * (dg[dSM3::uwv] + dg[dSM3::vwu] - dg[dSM3::uvw]),  // wuv
      0.5 * (dg[dSM3::uww] + dg[dSM3::wwu] - dg[dSM3::uww]),  // wuw
      0.5 * (dg[dSM3::vwv] + dg[dSM3::vwv] - dg[dSM3::vvw]),  // wvv
      0.5 * (dg[dSM3::vww] + dg[dSM3::wwv] - dg[dSM3::vww]),  // wvw
      0.5 * (dg[dSM3::www] + dg[dSM3::www] - dg[dSM3::www])  // www
  };
}

//! Inertial force @f$ F^k = - \Gamma^k_{ij} \, \dot{q}^i \, \dot{q}^j @f$.
inline IR3 metric_covariant::inertial_force(
    const IR3& q, const IR3& dot_q) const {
  ddIR3 gamma = christoffel_second_kind(q);
  SM3 dq2 = {
      dot_q[IR3::u] * dot_q[IR3::u], dot_q[IR3::u] * dot_q[IR3::v],
      dot_q[IR3::u] * dot_q[IR3::w], dot_q[IR3::v] * dot_q[IR3::v],
      dot_q[IR3::v] * dot_q[IR3::w], dot_q[IR3::w] * dot_q[IR3::w]};
 return {
    -gamma[ddIR3::uuu] * dq2[SM3::uu] - 2 * gamma[ddIR3::uuv] * dq2[SM3::uv] -
    2 * gamma[ddIR3::uuw] * dq2[SM3::uw] - gamma[ddIR3::uvv] * dq2[SM3::vv] -
    2 * gamma[ddIR3::uvw] * dq2[SM3::vw] - gamma[ddIR3::uww] * dq2[SM3::ww],
    -gamma[ddIR3::vuu] * dq2[SM3::uu] - 2 * gamma[ddIR3::vuv] * dq2[SM3::uv] -
    2 * gamma[ddIR3::vuw] * dq2[SM3::uw] - gamma[ddIR3::vvv] * dq2[SM3::vv] -
    2 * gamma[ddIR3::vvw] * dq2[SM3::vw] - gamma[ddIR3::vww] * dq2[SM3::ww],
    -gamma[ddIR3::wuu] * dq2[SM3::uu] - 2 * gamma[ddIR3::wuv] * dq2[SM3::uv] -
    2 * gamma[ddIR3::wuw] * dq2[SM3::uw] - gamma[ddIR3::wvv] * dq2[SM3::vv] -
    2 * gamma[ddIR3::wvw] * dq2[SM3::vw] - gamma[ddIR3::www] * dq2[SM3::ww]
 };
}

} // end namespace gyronimo.
