// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Manuel Assunção and Paulo Rodrigues.

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

// @metric_connected.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_connected.hh>

namespace gyronimo {

//! General covariant metric @f$g_{ij} = \mathbf{e}_i\cdot\mathbf{e}_j@f$.
SM3 metric_connected::operator()(const IR3& q) const {
  dIR3 e = my_morphism_->del(q);
  IR3 e1 = {e[dIR3::uu], e[dIR3::vu], e[dIR3::wu]};
  IR3 e2 = {e[dIR3::uv], e[dIR3::vv], e[dIR3::wv]};
  IR3 e3 = {e[dIR3::uw], e[dIR3::vw], e[dIR3::ww]};
  SM3 g = {inner_product(e1, e1), inner_product(e1, e2), inner_product(e1, e3),
           inner_product(e2, e2), inner_product(e2, e3), inner_product(e3, e3)};
  return g;
}

//! General covariant-metric derivatives from parent `morphism`.
/*!
    Extracts the metric derivatives from the Christoffel symbol of the first
    kind following the rule @f$\partial_k \, g_{ij} = \Gamma_{ijk} +
    \Gamma_{jik}@f$.
*/
dSM3 metric_connected::del(const IR3& q) const {
  ddIR3 gamma = christoffel_first_kind(q);
  return {
      gamma[ddIR3::uuu] + gamma[ddIR3::uuu],  // uuu
      gamma[ddIR3::uuv] + gamma[ddIR3::uuv],  // uuv
      gamma[ddIR3::uuw] + gamma[ddIR3::uuw],  // uuw
      gamma[ddIR3::uuv] + gamma[ddIR3::vuu],  // uvu
      gamma[ddIR3::uvv] + gamma[ddIR3::vuv],  // uvv
      gamma[ddIR3::uvw] + gamma[ddIR3::vuw],  // uvw
      gamma[ddIR3::uuw] + gamma[ddIR3::wuu],  // uwu
      gamma[ddIR3::uvw] + gamma[ddIR3::wuv],  // uwv
      gamma[ddIR3::uww] + gamma[ddIR3::wuw],  // uww
      gamma[ddIR3::vuv] + gamma[ddIR3::vuv],  // vvu
      gamma[ddIR3::vvv] + gamma[ddIR3::vvv],  // vvv
      gamma[ddIR3::vvw] + gamma[ddIR3::vvw],  // vvw
      gamma[ddIR3::vuw] + gamma[ddIR3::wuv],  // vwu
      gamma[ddIR3::vvw] + gamma[ddIR3::wvv],  // vwv
      gamma[ddIR3::vww] + gamma[ddIR3::wvw],  // vww
      gamma[ddIR3::wuw] + gamma[ddIR3::wuw],  // wwu
      gamma[ddIR3::wvw] + gamma[ddIR3::wvw],  // wwv
      gamma[ddIR3::www] + gamma[ddIR3::www]  // www
  };
}

//! General jacobian gradient from parent `morphism`.
/*!
    @f{equation*}{
    \partial_i J = J \left(
        \Gamma^1_{i 1} + \Gamma^2_{i 2} + \Gamma^3_{i 3} \right)
    @f}
*/
IR3 metric_connected::del_jacobian(const IR3& q) const {
  dIR3 ee = my_morphism_->del_inverse(q);
  ddIR3 de = my_morphism_->ddel(q);
  IR3 con = {
      ee[dIR3::uu] * de[ddIR3::uuu] + ee[dIR3::uv] * de[ddIR3::vuu] +
          ee[dIR3::uw] * de[ddIR3::wuu] + ee[dIR3::vu] * de[ddIR3::uuv] +
          ee[dIR3::vv] * de[ddIR3::vuv] + ee[dIR3::vw] * de[ddIR3::wuv] +
          ee[dIR3::wu] * de[ddIR3::uuw] + ee[dIR3::wv] * de[ddIR3::vuw] +
          ee[dIR3::ww] * de[ddIR3::wuw],
      ee[dIR3::uu] * de[ddIR3::uuv] + ee[dIR3::uv] * de[ddIR3::vuv] +
          ee[dIR3::uw] * de[ddIR3::wuv] + ee[dIR3::vu] * de[ddIR3::uvv] +
          ee[dIR3::vv] * de[ddIR3::vvv] + ee[dIR3::vw] * de[ddIR3::wvv] +
          ee[dIR3::wu] * de[ddIR3::uvw] + ee[dIR3::wv] * de[ddIR3::vvw] +
          ee[dIR3::ww] * de[ddIR3::wvw],
      ee[dIR3::uu] * de[ddIR3::uuw] + ee[dIR3::uv] * de[ddIR3::vuw] +
          ee[dIR3::uw] * de[ddIR3::wuw] + ee[dIR3::vu] * de[ddIR3::uvw] +
          ee[dIR3::vv] * de[ddIR3::vvw] + ee[dIR3::vw] * de[ddIR3::wvw] +
          ee[dIR3::wu] * de[ddIR3::uww] + ee[dIR3::wv] * de[ddIR3::vww] +
          ee[dIR3::ww] * de[ddIR3::www]};
  return (jacobian(q) * con);
}

}  // end namespace gyronimo.
