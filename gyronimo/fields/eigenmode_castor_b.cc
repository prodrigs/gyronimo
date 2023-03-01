// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues.

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

// @eigenmode_castor_b.cc, this file is part of ::gyronimo::

#include <ranges>
#include <gyronimo/fields/eigenmode_castor_b.hh>

namespace gyronimo{

eigenmode_castor_b::eigenmode_castor_b(
    double m_factor, double t_factor,
    const parser_castor *p, const metric_helena *g,
    const interpolator1d_factory* ifactory)
    : IR3field_c1(m_factor, t_factor, g),
      norm_factor_(1.0),
      parser_(p), metric_(g),
      w_(p->eigenvalue_real(), p->eigenvalue_imag()),
      n_tor_squared_(p->n_tor()*p->n_tor()),
      i_n_tor_(0.0, p->n_tor()),
      tildeA1_(p->s(), p->a1_real(), p->a1_imag(), p->m(), ifactory),
      tildeA2_(p->s(), p->a2_real(), p->a2_imag(), p->m(), ifactory),
      tildeA3_(p->s(), p->a3_real(), p->a3_imag(), p->m(), ifactory) {
  using namespace std;
  auto max_magnitude_at_radius = [this](double s) {
    auto highest_harmonic = ranges::max(views::transform(
        this->parser_->m(), [](auto m) {return abs(m);}));
    size_t chi_range_size = (size_t)(8*highest_harmonic);
    double delta_chi = 2*numbers::pi/chi_range_size;
    auto chi_range = views::transform(
        views::iota(0u, chi_range_size), bind_front(multiplies{}, delta_chi));
    return ranges::max(views::transform(
        chi_range,
        [this, s](double chi){return this->magnitude({s, chi, 0}, 0);}));
  };
  norm_factor_ = 1.0/ranges::max(
      views::transform(parser_->s(), max_magnitude_at_radius));
}

//! Magnetic field \f$B^i = \epsilon^{ijk}/\sqrt{g} \partial_j A_k\f$.
IR3 eigenmode_castor_b::contravariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double phi = position[IR3::w];
  double chi = metric_->reduce_chi(position[IR3::v]);
  std::complex<double> factor = norm_factor_*
      this->exp_wt_nphi(time, phi)/this->metric()->jacobian(position);
  return {
      std::real(factor*(this->d2A3(s, chi) - this->d3A2(s, chi))),
      std::real(factor*(this->d3A1(s, chi) - this->d1A3(s, chi))),
      std::real(factor*(this->d1A2(s, chi) - this->d2A1(s, chi)))};
}

//! Magnetic-field spatial derivatives.
/*!
    Implements the formula
    \f[ \partial_l B^i = \frac{1}{\sqrt{g}} \bigl[
        \epsilon^{ijk} \partial^2_{jl} A_k - B^i \partial_l \sqrt{g} \bigr] \f]
*/
dIR3 eigenmode_castor_b::del_contravariant(
    const IR3& position, double time) const {
  double s = position[IR3::u];
  double chi = metric_->reduce_chi(position[IR3::v]);
  double phi = position[IR3::w];
  std::complex<double> factor = norm_factor_*this->exp_wt_nphi(time, phi);
  std::valarray epsilon_ijk_partial2_jl_A_k = {
      std::real(factor*(this->d21A3(s, chi) - this->d31A2(s, chi))),
      std::real(factor*(this->d22A3(s, chi) - this->d32A2(s, chi))),
      std::real(factor*(this->d23A3(s, chi) - this->d33A2(s, chi))),
      std::real(factor*(this->d31A1(s, chi) - this->d11A3(s, chi))),
      std::real(factor*(this->d32A1(s, chi) - this->d12A3(s, chi))),
      std::real(factor*(this->d33A1(s, chi) - this->d13A3(s, chi))),
      std::real(factor*(this->d11A2(s, chi) - this->d21A1(s, chi))),
      std::real(factor*(this->d12A2(s, chi) - this->d22A1(s, chi))),
      std::real(factor*(this->d13A2(s, chi) - this->d23A1(s, chi)))};
  IR3 B = this->contravariant(position, time);
  IR3 dg = this->metric()->del_jacobian(position);
  std::valarray B_i_partial_l_sqrt_g = {
      B[IR3::u]*dg[IR3::u],B[IR3::u]*dg[IR3::v],B[IR3::u]*dg[IR3::w],
      B[IR3::v]*dg[IR3::u],B[IR3::v]*dg[IR3::v],B[IR3::v]*dg[IR3::w],
      B[IR3::w]*dg[IR3::u],B[IR3::w]*dg[IR3::v],B[IR3::w]*dg[IR3::w]};
  double ijacobian = 1.0/this->metric()->jacobian(position);
  dIR3 result;
  result = ijacobian*(epsilon_ijk_partial2_jl_A_k - B_i_partial_l_sqrt_g);
  return result;
}

//! Time derivative \f$\partial_t\mathbf{B} = i\omega\nabla\times\mathbf{A}\f$.
IR3 eigenmode_castor_b::partial_t_contravariant(
    const IR3& position, double time) const {
  double s = position[IR3::u];
  double phi = position[IR3::w];
  double chi = metric_->reduce_chi(position[IR3::v]);
  std::complex<double> factor = norm_factor_*this->exp_wt_nphi(time, phi);
  using namespace std::complex_literals;
  factor *= 1i*w_/this->metric()->jacobian(position);
  return {
      std::real(factor*(this->d2A3(s, chi) - this->d3A2(s, chi))),
      std::real(factor*(this->d3A1(s, chi) - this->d1A3(s, chi))),
      std::real(factor*(this->d1A2(s, chi) - this->d2A1(s, chi)))};
}

inline
std::complex<double> eigenmode_castor_b::exp_wt_nphi(
    double time, double phi) const {
  using namespace std::complex_literals;
  return std::exp(w_*time + 1i*parser_->n_tor()*phi);
}
inline
std::complex<double> eigenmode_castor_b::d1A2(double s, double chi) const {
  return tildeA2_.partial_u(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d1A3(double s, double chi) const {
  return tildeA3_.partial_u(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d2A1(double s, double chi) const {
  return tildeA1_.partial_v(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d2A3(double s, double chi) const {
  return tildeA3_.partial_v(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d3A1(double s, double chi) const {
  using namespace std::complex_literals;
  return i_n_tor_*tildeA1_(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d3A2(double s, double chi) const {
  using namespace std::complex_literals;
  return i_n_tor_*tildeA2_(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d11A2(double s, double chi) const {
  return tildeA2_.partial2_uu(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d11A3(double s, double chi) const {
  return tildeA3_.partial2_uu(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d21A1(double s, double chi) const {
  return tildeA1_.partial2_uv(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d21A3(double s, double chi) const {
  return tildeA3_.partial2_uv(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d31A1(double s, double chi) const {
  return i_n_tor_*tildeA1_.partial_u(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d31A2(double s, double chi) const {
  return i_n_tor_*tildeA2_.partial_u(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d12A2(double s, double chi) const {
  return tildeA2_.partial2_uv(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d12A3(double s, double chi) const {
  return tildeA3_.partial2_uv(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d22A1(double s, double chi) const {
  return tildeA1_.partial2_vv(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d22A3(double s, double chi) const {
  return tildeA3_.partial2_vv(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d32A1(double s, double chi) const {
  return i_n_tor_*tildeA1_.partial_v(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d32A2(double s, double chi) const {
  return i_n_tor_*tildeA2_.partial_v(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d13A2(double s, double chi) const {
  return i_n_tor_*tildeA2_.partial_u(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d13A3(double s, double chi) const {
  return i_n_tor_*tildeA3_.partial_u(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d23A1(double s, double chi) const {
  return i_n_tor_*tildeA1_.partial_v(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d23A3(double s, double chi) const {
  return i_n_tor_*tildeA3_.partial_v(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d33A1(double s, double chi) const {
  return -n_tor_squared_*tildeA1_(s, chi);
}
inline
std::complex<double> eigenmode_castor_b::d33A2(double s, double chi) const {
  return -n_tor_squared_*tildeA2_(s, chi);
}

} // end namespace gyronimo.
