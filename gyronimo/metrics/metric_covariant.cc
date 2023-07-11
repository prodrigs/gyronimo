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

#include <gyronimo/core/contraction.hh>
#include <gyronimo/metrics/metric_covariant.hh>

#include <cmath>

namespace gyronimo {

double metric_covariant::jacobian(const IR3& r) const {
  SM3 g = (*this)(r);
  return std::sqrt(
      g[SM3::uu] * g[SM3::vv] * g[SM3::ww] +
      2.0 * g[SM3::uv] * g[SM3::uw] * g[SM3::vw] -
      g[SM3::uv] * g[SM3::uv] * g[SM3::ww] -
      g[SM3::uu] * g[SM3::vw] * g[SM3::vw] -
      g[SM3::uw] * g[SM3::uw] * g[SM3::vv]);
}
IR3 metric_covariant::del_jacobian(const IR3& r) const {
  SM3 g = (*this)(r);
  dSM3 dg = this->del(r);
  double ijacobian = 1.0 / this->jacobian(r);

  IR3 aux_1a = {
      g[SM3::uw] * g[SM3::vw], g[SM3::uv] * g[SM3::vw],
      g[SM3::uv] * g[SM3::uw]};
  dIR3 aux_1b = {dg[dSM3::uvu], dg[dSM3::uwu], dg[dSM3::vwu],
                 dg[dSM3::uvv], dg[dSM3::uwv], dg[dSM3::vwv],
                 dg[dSM3::uvw], dg[dSM3::uww], dg[dSM3::vww]};
  IR3 aux_1 = contraction<second>(aux_1b, aux_1a);

  IR3 aux_2a = {
      g[SM3::vv] * g[SM3::ww], g[SM3::uu] * g[SM3::ww],
      g[SM3::uu] * g[SM3::vv]};
  dIR3 aux_2b = {dg[dSM3::uuu], dg[dSM3::vvu], dg[dSM3::wwu],
                 dg[dSM3::uuv], dg[dSM3::vvv], dg[dSM3::wwv],
                 dg[dSM3::uuw], dg[dSM3::vvw], dg[dSM3::www]};
  IR3 aux_2 = contraction<second>(aux_2b, aux_2a);

  IR3 aux_3a = {
      g[SM3::uu] * g[SM3::vw], g[SM3::vv] * g[SM3::uw],
      g[SM3::ww] * g[SM3::uv]};
  dIR3 aux_3b = {dg[dSM3::vwu], dg[dSM3::uwu], dg[dSM3::uvu],
                 dg[dSM3::vwv], dg[dSM3::uwv], dg[dSM3::uvv],
                 dg[dSM3::vww], dg[dSM3::uww], dg[dSM3::uvw]};
  IR3 aux_3 = contraction<second>(aux_3b, aux_3a);

  IR3 aux_4a = {
      g[SM3::vw] * g[SM3::vw], g[SM3::uw] * g[SM3::uw],
      g[SM3::uv] * g[SM3::uv]};
  IR3 aux_4 = contraction<second>(aux_2b, aux_4a);

  return ijacobian * (aux_1 + 0.5 * aux_2 - aux_3 - 0.5 * aux_4);
}
SM3 metric_covariant::inverse(const IR3& r) const {
  SM3 m = (*this)(r);
  return gyronimo::inverse(m);
}

//! General-purpose product of a covariant metric and a contravariant vector.
/*!
    Returns the *covariant* product of the *covariant* metric, evaluated at a
    contravariant position, by the *contravariant* vector `B`.
*/
IR3 metric_covariant::to_covariant(const IR3& B, const IR3& r) const {
  return contraction((*this)(r), B);
}

//! General-purpose product of a contravariant metric and a covariant vector.
/*!
    Returns the *contravariant* product of the *contravariant* metric, evaluated
    at a contravariant position, with the *covariant* vector `B`. This method is
    significantly more expensive than `to_covariant` because it involves the
    covariant-metric inversion.
*/
IR3 metric_covariant::to_contravariant(const IR3& B, const IR3& r) const {
  return contraction(this->inverse(r), B);
}

//! General-purpose derivative of the inverse metric.
/*!
    Implements the rule @f$
    \partial_k g^{ij} = - g^{im} \partial_k g_{mn} g^{nj}
    @f$.
*/
dSM3 metric_covariant::del_inverse(const IR3& r) const {
  SM3 ig = this->inverse(r);
  dSM3 dg = contraction(ig, this->del(r), ig);
  return {-dg[dSM3::uuu], -dg[dSM3::uuv], -dg[dSM3::uuw], -dg[dSM3::uvu],
          -dg[dSM3::uvv], -dg[dSM3::uvw], -dg[dSM3::uwu], -dg[dSM3::uwv],
          -dg[dSM3::uww], -dg[dSM3::vvu], -dg[dSM3::vvv], -dg[dSM3::vvw],
          -dg[dSM3::vwu], -dg[dSM3::vwv], -dg[dSM3::vww], -dg[dSM3::wwu],
          -dg[dSM3::wwv], -dg[dSM3::www]};
}

//! General-purpose Christoffel symbols of the first kind.
/*!
    Implements the rule @f$
    \Gamma_{ijk} = \frac{1}{2} \left(
                \frac{\partial g_{ij}}{\partial q^k} +
                \frac{\partial g_{ik}}{\partial q^j} -
                \frac{\partial g_{jk}}{\partial q^i} \right)
    @f$.
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

//! General-purpose Christoffel symbols of the second kind.
/*!
    Implements the rule @f$ \Gamma^k_{ij} = g^{km} \, \Gamma_{mij} @f$.
*/
ddIR3 metric_covariant::christoffel_second_kind(const IR3& q) const {
  return contraction<second, first>(
      this->inverse(q), this->christoffel_first_kind(q));
}

//! Contravariant inertial force.
/*!
    Implements the rule @f$ F^k = - \Gamma^k_{ij} \, \dot{q}^i \, \dot{q}^j @f$.
*/
IR3 metric_covariant::inertial_force(const IR3& q, const IR3& dot_q) const {
  ddIR3 gamma = christoffel_second_kind(q);
  return (-1.0)*contraction<second,third>(gamma, dot_q, dot_q);
}

} // end namespace gyronimo.
