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

// @eigenmode_castor_a.cc, this file is part of ::gyronimo::

#include <ranges>
#include <gyronimo/fields/eigenmode_castor_a.hh>

namespace gyronimo{

eigenmode_castor_a::eigenmode_castor_a(
    double m_factor, double t_factor,
    const parser_castor *p, const metric_helena *g,
    const interpolator1d_factory* ifactory)
    : IR3field(m_factor, t_factor, g),
      norm_factor_(1.0), parser_(p), metric_(g),
      w_(p->eigenvalue_real(), p->eigenvalue_imag()), i_n_tor_(0.0, p->n_tor()),
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

eigenmode_castor_a::eigenmode_castor_a(
    double m_factor, double t_factor,
    const parser_castor *p, const metric_helena *g,
    const interpolator1d_factory* ifactory, double norm_factor)
    : IR3field(m_factor, t_factor, g),
      norm_factor_(norm_factor), parser_(p), metric_(g),
      w_(p->eigenvalue_real(), p->eigenvalue_imag()), i_n_tor_(0.0, p->n_tor()),
      tildeA1_(p->s(), p->a1_real(), p->a1_imag(), p->m(), ifactory),
      tildeA2_(p->s(), p->a2_real(), p->a2_imag(), p->m(), ifactory),
      tildeA3_(p->s(), p->a3_real(), p->a3_imag(), p->m(), ifactory) {
}

IR3 eigenmode_castor_a::covariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double phi = position[IR3::w];
  double chi = metric_->reduce_chi(position[IR3::v]);
  std::complex<double> factor = norm_factor_*std::exp(w_*time + i_n_tor_*phi);
  return {
      std::real(factor*(this->tildeA1_(s, chi))),
          std::real(factor*(this->tildeA2_(s, chi))),
              std::real(factor*(this->tildeA3_(s, chi)))};
}
IR3 eigenmode_castor_a::contravariant(const IR3& position, double time) const {
  return this->metric()->
      to_contravariant(this->covariant(position, time), position);
}

} // end namespace gyronimo.
