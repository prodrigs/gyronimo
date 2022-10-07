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

// @metric_covariant.cc, this file is part of ::gyronimo::

#include <cmath>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! General-purpose implementation of the Jacobian.
double metric_covariant::jacobian(const IR3& r) const {
  SM3 g = (*this)(r);
  return std::sqrt(g[SM3::uu]*g[SM3::vv]*g[SM3::ww] +
      2.0*g[SM3::uv]*g[SM3::uw]*g[SM3::vw] - g[SM3::uv]*g[SM3::uv]*g[SM3::ww] -
          g[SM3::uu]*g[SM3::vw]*g[SM3::vw] - g[SM3::uw]*g[SM3::uw]*g[SM3::vv]);
}

//! General-purpose implementation of the Jacobian gradient.
IR3 metric_covariant::del_jacobian(const IR3& r) const {
  SM3 g = (*this)(r);
  dSM3 dg = this->del(r);
  return {
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
}

//! General-purpose product of a covariant metric and a contravariant vector.
/*!
    Returns the *covariant* product of the *covariant* metric, evaluated at a
    contravariant position, by the *contravariant* vector `B`.
*/
IR3 metric_covariant::to_covariant(const IR3& B, const IR3& r) const {
  return contraction((*this)(r), B);
}

//! General-purpose implementation of the inverse (i.e., contravariant metric).
SM3 metric_covariant::inverse(const IR3& r) const {
  SM3 m = (*this)(r);
  return gyronimo::inverse(m);
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

//! General-purpose implementation of the inverse derivatives.
/*!
    Implements the rule
    @f$ \partial_k g^{ij} = - g^{im} \partial_k g_{mn} g^{nj} @f$
*/
dSM3 metric_covariant::del_inverse(const IR3& r) const {
  SM3 ig = this->inverse(r);
  dSM3 dg = contraction(ig, this->del(r), ig);
  return {
    -dg[dSM3::uuu], -dg[dSM3::uuv], -dg[dSM3::uuw],
    -dg[dSM3::uvu], -dg[dSM3::uvv], -dg[dSM3::uvw],
    -dg[dSM3::uwu], -dg[dSM3::uwv], -dg[dSM3::uww],
    -dg[dSM3::vvu], -dg[dSM3::vvv], -dg[dSM3::vvw],
    -dg[dSM3::vwu], -dg[dSM3::vwv], -dg[dSM3::vww],
    -dg[dSM3::wwu], -dg[dSM3::wwv], -dg[dSM3::www]};
}

//! Implementation of Christoffel symbols of the first kind @f$ \Gamma_{\alpha\beta\gamma} @f$
/*!
    Implements the rule
    @f$ \Gamma_{\alpha\beta\gamma} = \frac{1}{2} \left( 
		\frac{\partial g_{\alpha\beta}}{\partial q^\gamma} +
		\frac{\partial g_{\alpha\gamma}}{\partial q^\beta} -
		\frac{\partial g_{\beta\gamma}}{\partial q^\alpha}
	\right) @f$
*/
ddIR3 metric_covariant::christoffel_first_kind(const IR3& r) const {
	dSM3 dg = del(r);
	return {
		0.5 * (dg[dSM3::uuu] + dg[dSM3::uuu] - dg[dSM3::uuu]), // uuu
		0.5 * (dg[dSM3::uuv] + dg[dSM3::uvu] - dg[dSM3::uvu]), // uuv
		0.5 * (dg[dSM3::uuw] + dg[dSM3::uwu] - dg[dSM3::uwu]), // uuw
		0.5 * (dg[dSM3::uvv] + dg[dSM3::uvv] - dg[dSM3::vvu]), // uvv
		0.5 * (dg[dSM3::uvw] + dg[dSM3::uwv] - dg[dSM3::vwu]), // uvw
		0.5 * (dg[dSM3::uww] + dg[dSM3::uww] - dg[dSM3::wwu]), // uww
		0.5 * (dg[dSM3::uvu] + dg[dSM3::uvu] - dg[dSM3::uuv]), // vuu
		0.5 * (dg[dSM3::uvv] + dg[dSM3::vvu] - dg[dSM3::uvv]), // vuv
		0.5 * (dg[dSM3::uvw] + dg[dSM3::vwu] - dg[dSM3::uwv]), // vuw
		0.5 * (dg[dSM3::vvv] + dg[dSM3::vvv] - dg[dSM3::vvv]), // vvv
		0.5 * (dg[dSM3::vvw] + dg[dSM3::vwv] - dg[dSM3::vwv]), // vvw
		0.5 * (dg[dSM3::vww] + dg[dSM3::vww] - dg[dSM3::wwv]), // vww
    	0.5 * (dg[dSM3::uwu] + dg[dSM3::uwu] - dg[dSM3::uuw]), // wuu
		0.5 * (dg[dSM3::uwv] + dg[dSM3::vwu] - dg[dSM3::uvw]), // wuv
		0.5 * (dg[dSM3::uww] + dg[dSM3::wwu] - dg[dSM3::uww]), // wuw
		0.5 * (dg[dSM3::vwv] + dg[dSM3::vwv] - dg[dSM3::vvw]), // wvv
		0.5 * (dg[dSM3::vww] + dg[dSM3::wwv] - dg[dSM3::vww]), // wvw
		0.5 * (dg[dSM3::www] + dg[dSM3::www] - dg[dSM3::www])  // www
	};
}

//! Implementation of Christoffel symbols of the second kind @f$ \Gamma^\alpha_{\beta\gamma} @f$
/*!
    Implements the rule
    @f$ \Gamma^\alpha_{\beta\gamma} = g^{\alpha\mu} \, \Gamma_{\mu\beta\gamma} @f$
*/
ddIR3 metric_covariant::christoffel_second_kind(const IR3& r) const {
	SM3 ig = inverse(r);
	ddIR3 CF1 = christoffel_first_kind(r);
	return contraction<second, first>(ig, CF1);
}

} // end namespace gyronimo.
