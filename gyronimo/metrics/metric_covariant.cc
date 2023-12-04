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
      2.0*g[SM3::uv]*g[SM3::uw]*g[SM3::vw] - g[SM3::uv]*g[SM3::uv]*g[SM3::ww] -
      g[SM3::uu]*g[SM3::vw]*g[SM3::vw] - g[SM3::uw]*g[SM3::uw]*g[SM3::vv];
  return std::sqrt(det_g);
}

//! General-purpose implementation of the Jacobian gradient.
IR3 metric_covariant::del_jacobian(const IR3& r) const {
  SM3 g = (*this)(r);
  dSM3 dg = this->del(r);
  double det_g = g[SM3::uu]*g[SM3::vv]*g[SM3::ww] +
      2.0*g[SM3::uv]*g[SM3::uw]*g[SM3::vw] - g[SM3::uv]*g[SM3::uv]*g[SM3::ww] -
      g[SM3::uu]*g[SM3::vw]*g[SM3::vw] - g[SM3::uw]*g[SM3::uw]*g[SM3::vv];
  IR3 grad_det_g = {
      2.0*g[SM3::vw]*(
          g[SM3::uw]*dg[dSM3::uvu] +
          g[SM3::uv]*dg[dSM3::uwu] - g[SM3::uu]*dg[dSM3::vwu]) -
      g[SM3::uv]*g[SM3::uv]*dg[dSM3::wwu] + g[SM3::vv]*(
          g[SM3::ww]*dg[dSM3::uuu] -
          2.0*g[SM3::uw]*dg[dSM3::uwu] + g[SM3::uu]*dg[dSM3::wwu]) -
      g[SM3::vw]*g[SM3::vw]*dg[dSM3::uuu] -
      2.0*g[SM3::uv]*g[SM3::ww]*dg[dSM3::uvu] - 
      g[SM3::uw]*g[SM3::uw]*dg[dSM3::vvu] +
      g[SM3::uu]*g[SM3::ww]*dg[dSM3::vvu] + 
      2.0*g[SM3::uv]*g[SM3::uw]*dg[dSM3::vwu],
      2.0*g[SM3::vw]*(
          g[SM3::uw]*dg[dSM3::uvv] +
          g[SM3::uv]*dg[dSM3::uwv] - g[SM3::uu]*dg[dSM3::vwv]) -
      g[SM3::uv]*g[SM3::uv]*dg[dSM3::wwv] + g[SM3::vv]*(
          g[SM3::ww]*dg[dSM3::uuv] -
          2.0*g[SM3::uw]*dg[dSM3::uwv] + g[SM3::uu]*dg[dSM3::wwv]) -
      g[SM3::vw]*g[SM3::vw]*dg[dSM3::uuv] -
      2.0*g[SM3::uv]*g[SM3::ww]*dg[dSM3::uvv] - 
      g[SM3::uw]*g[SM3::uw]*dg[dSM3::vvv] +
      g[SM3::uu]*g[SM3::ww]*dg[dSM3::vvv] + 
      2.0*g[SM3::uv]*g[SM3::uw]*dg[dSM3::vwv],
      2.0*g[SM3::vw]*(
          g[SM3::uw]*dg[dSM3::uvw] +
          g[SM3::uv]*dg[dSM3::uww] - g[SM3::uu]*dg[dSM3::vww]) -
      g[SM3::uv]*g[SM3::uv]*dg[dSM3::www] + g[SM3::vv]*(
          g[SM3::ww]*dg[dSM3::uuw] -
          2.0*g[SM3::uw]*dg[dSM3::uww] + g[SM3::uu]*dg[dSM3::www]) -
      g[SM3::vw]*g[SM3::vw]*dg[dSM3::uuw] -
      2.0*g[SM3::uv]*g[SM3::ww]*dg[dSM3::uvw] - 
      g[SM3::uw]*g[SM3::uw]*dg[dSM3::vvw] +
      g[SM3::uu]*g[SM3::ww]*dg[dSM3::vvw] + 
      2.0*g[SM3::uv]*g[SM3::uw]*dg[dSM3::vww]};
  return (0.5/std::sqrt(det_g))*grad_det_g;
}

//! General derivative @f$\partial_k g^{ij} = -g^{im}\partial_k g_{mn}g^{nj}@f$.
dSM3 metric_covariant::del_inverse(const IR3& q) const {
  SM3 ig = this->inverse(q);
  dSM3 dg = contraction(ig, this->del(q), ig);
  return {-dg[dSM3::uuu], -dg[dSM3::uuv], -dg[dSM3::uuw], -dg[dSM3::uvu],
          -dg[dSM3::uvv], -dg[dSM3::uvw], -dg[dSM3::uwu], -dg[dSM3::uwv],
          -dg[dSM3::uww], -dg[dSM3::vvu], -dg[dSM3::vvv], -dg[dSM3::vvw],
          -dg[dSM3::vwu], -dg[dSM3::vwv], -dg[dSM3::vww], -dg[dSM3::wwu],
          -dg[dSM3::wwv], -dg[dSM3::www]};
}

//! General Christoffel symbol @f$\Gamma_{ijk}@f$.
/*! @f{equation*}{
      \Gamma_{ijk} = \frac{1}{2} \left(
                \frac{\partial g_{ij}}{\partial q^k} +
                \frac{\partial g_{ik}}{\partial q^j} -
                \frac{\partial g_{jk}}{\partial q^i} \right)
    @f}
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

} // end namespace gyronimo.
