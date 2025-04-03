// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @IR3field_c1.cc, this file is part of ::gyronimo::

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/core/contraction.hh>

namespace gyronimo {

//! Partial derivatives of the covariant components.
/*!
    Implements the rule
    \f$\partial_i E_j = \partial_i(g_{jk} E^k)
      = g_{jk} \partial_i(E^k) + \partial_i(g_{jk}) E^k\f$
*/
dIR3 IR3field_c1::del_covariant(const IR3& position, double time) const {
  const metric_covariant *g = this->metric();
  dIR3 c1 = contraction<second>(
      (*g).del(position), this->contravariant(position, time));
  dIR3 c2 = contraction<first>(
      this->del_contravariant(position, time), (*g)(position));
  return {
      c1[dIR3::uu] + c2[dIR3::uu], c1[dIR3::uv] + c2[dIR3::uv],
      c1[dIR3::uw] + c2[dIR3::uw], c1[dIR3::vu] + c2[dIR3::vu],
      c1[dIR3::vv] + c2[dIR3::vv], c1[dIR3::vw] + c2[dIR3::vw],
      c1[dIR3::wu] + c2[dIR3::wu], c1[dIR3::wv] + c2[dIR3::wv],
      c1[dIR3::ww] + c2[dIR3::ww]};
}

//! Partial time derivatives of the covariant components.
/*!
    Implements the rule
    \f$\partial_t E_j = \partial_t(g_{jk} E^k) = g_{jk} \partial_t(E^k)\f$
*/
IR3 IR3field_c1::partial_t_covariant(const IR3& position, double time) const {
  IR3 dE = this->partial_t_contravariant(position, time);
  return this->metric()->to_covariant(dE, position);
}

//! Divergence of the field.
/*!
    Implements the rule
    \f$ \nabla\cdot\mathbf{E} =
        \partial_i E^i + \sqrt{g}^{-1} E^i \partial_i \sqrt{g} \f$
*/
double IR3field_c1::div(const IR3& position, double time) const {
  double ijacobian = 1.0 / this->metric()->jacobian(position);
  IR3 E = this->contravariant(position, time);
  dIR3 dE = this->del_contravariant(position, time);
  return dE[dIR3::uu] + dE[dIR3::vv] + dE[dIR3::ww] +
      ijacobian * inner_product(E, this->metric()->del_jacobian(position));
}

//! Contravariant components of the curl operator.
/*!
    Implements the rule
    \f$ \sqrt{g}(\nabla\times\mathbf{E})^k = \epsilon^{kij}\partial_iE_j \f$
*/
IR3 IR3field_c1::curl(const IR3& position, double time) const {
  double ijacobian = 1.0 / this->metric()->jacobian(position);
  dIR3 dE = this->del_covariant(position, time);
  return {
      (dE[dIR3::wv] - dE[dIR3::vw])*ijacobian,
      (dE[dIR3::uw] - dE[dIR3::wu])*ijacobian,
      (dE[dIR3::vu] - dE[dIR3::uv])*ijacobian};
}

//! Covariant components of the magnitude gradient.
/*!
    Implements the rules
    \f$B^2 = B_j B^j; \quad
       2 B \partial_i B = 2 B (B_j \partial_i B^j + B^j \partial_i B_j)\f$
*/
IR3 IR3field_c1::del_magnitude(const IR3& position, double time) const {
  return (0.5/this->magnitude(position, time))*(
      contraction<first>(
          this->del_covariant(position, time),
          this->contravariant(position, time)) +
      contraction<first>(
          this->del_contravariant(position, time),
          this->covariant(position, time)));
}

//! Partial time derivative of the magnitude.
/*!
    Implements the rules
    \f$B^2 = B_j B^j; \quad
       2 B \partial_t B = 2 B (B_j \partial_t B^j + B^j \partial_t B_j)\f$
*/
double IR3field_c1::partial_t_magnitude(
    const IR3& position, double time) const {
  return (0.5/this->magnitude(position, time))*(
      inner_product(
          this->partial_t_covariant(position, time),
          this->contravariant(position, time)) +
      inner_product(
          this->partial_t_contravariant(position, time),
          this->covariant(position, time)));
}

dIR3 versor_field::del_contravariant(const IR3& q, double t) const {
  IR3 Bj = base_field_->contravariant(q, t);
  IR3 dB = base_field_->del_magnitude(q, t);
  double iB = 1 / base_field_->magnitude(q, t);
  dIR3 dBj = base_field_->del_contravariant(q, t);
  return {
    iB * (dBj[dIR3::uu] - iB * Bj[IR3::u] * dB[IR3::u]),
    iB * (dBj[dIR3::uv] - iB * Bj[IR3::u] * dB[IR3::v]),
    iB * (dBj[dIR3::uw] - iB * Bj[IR3::u] * dB[IR3::w]),
    iB * (dBj[dIR3::vu] - iB * Bj[IR3::v] * dB[IR3::u]),
    iB * (dBj[dIR3::vv] - iB * Bj[IR3::v] * dB[IR3::v]),
    iB * (dBj[dIR3::vw] - iB * Bj[IR3::v] * dB[IR3::w]),
    iB * (dBj[dIR3::wu] - iB * Bj[IR3::w] * dB[IR3::u]),
    iB * (dBj[dIR3::wv] - iB * Bj[IR3::w] * dB[IR3::v]),
    iB * (dBj[dIR3::ww] - iB * Bj[IR3::w] * dB[IR3::w])};
}

// Covariant derivate of the magnetic field versor
dIR3 versor_field::del_covariant(const IR3& q, double t) const {
  IR3 b_con = this->contravariant(q, t);
  dIR3 del_b_con = this->del_contravariant(q, t);
  ddIR3 christoffel = this->metric()->christoffel_second_kind(q);
  return {
    del_b_con[dIR3::uu] + christoffel[ddIR3::uuu]*b_con[IR3::u] + christoffel[ddIR3::uuv]*b_con[IR3::v] + christoffel[ddIR3::uuw]*b_con[IR3::w],
    del_b_con[dIR3::uv] + christoffel[ddIR3::uuv]*b_con[IR3::u] + christoffel[ddIR3::uvv]*b_con[IR3::v] + christoffel[ddIR3::uvw]*b_con[IR3::w],
    del_b_con[dIR3::uw] + christoffel[ddIR3::uuw]*b_con[IR3::u] + christoffel[ddIR3::uvw]*b_con[IR3::v] + christoffel[ddIR3::uww]*b_con[IR3::w],
    del_b_con[dIR3::vu] + christoffel[ddIR3::vuu]*b_con[IR3::u] + christoffel[ddIR3::vuv]*b_con[IR3::v] + christoffel[ddIR3::vuw]*b_con[IR3::w],
    del_b_con[dIR3::vv] + christoffel[ddIR3::vuv]*b_con[IR3::u] + christoffel[ddIR3::vvv]*b_con[IR3::v] + christoffel[ddIR3::vvw]*b_con[IR3::w],
    del_b_con[dIR3::vw] + christoffel[ddIR3::vuw]*b_con[IR3::u] + christoffel[ddIR3::vvw]*b_con[IR3::v] + christoffel[ddIR3::vww]*b_con[IR3::w],
    del_b_con[dIR3::wu] + christoffel[ddIR3::wuu]*b_con[IR3::u] + christoffel[ddIR3::wuv]*b_con[IR3::v] + christoffel[ddIR3::wuw]*b_con[IR3::w],
    del_b_con[dIR3::wv] + christoffel[ddIR3::wuv]*b_con[IR3::u] + christoffel[ddIR3::wvv]*b_con[IR3::v] + christoffel[ddIR3::wvw]*b_con[IR3::w],
    del_b_con[dIR3::ww] + christoffel[ddIR3::wuw]*b_con[IR3::u] + christoffel[ddIR3::wvw]*b_con[IR3::v] + christoffel[ddIR3::www]*b_con[IR3::w]};
}

} // end namespace gyronimo.
