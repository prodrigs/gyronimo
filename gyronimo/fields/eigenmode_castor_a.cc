// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues and João Palma.

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

//! Normalises native `castor` harmonics to the poloidal-section maximum.
eigenmode_castor_a::eigenmode_castor_a(
    double m_factor, double v_alfven,
    const parser_castor *p, const metric_helena *g,
    const interpolator1d_factory* ifactory)
    : IR3field(m_factor, g->parser()->rmag()/v_alfven, g),
      native_factor_(1.0), parser_(p), metric_(g),
      eigenvalue_(p->eigenvalue_real(), p->eigenvalue_imag()),
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
  native_factor_ = 1.0/ranges::max(
      views::transform(parser_->s(), max_magnitude_at_radius));
}

//! Builds with normalisation compatible with a parent `eigenmode_castor_b`.
eigenmode_castor_a::eigenmode_castor_a(
    const eigenmode_castor_b* parent, const interpolator1d_factory* ifactory)
    : IR3field(
          parent->m_factor()*(parent->v_alfven()*parent->t_factor()),
          parent->t_factor(),
          static_cast<const metric_helena*>(parent->metric())),
      native_factor_(parent->native_factor()),
      parser_(parent->parser()),
      metric_(static_cast<const metric_helena*>(parent->metric())),
      eigenvalue_(parent->parser()->eigenvalue_real(),
          parent->parser()->eigenvalue_imag()),
      i_n_tor_(0.0, parent->parser()->n_tor()),
      tildeA1_(parent->parser()->s(), parent->parser()->a1_real(),
          parent->parser()->a1_imag(), parent->parser()->m(), ifactory),
      tildeA2_(parent->parser()->s(), parent->parser()->a2_real(),
          parent->parser()->a2_imag(), parent->parser()->m(), ifactory),
      tildeA3_(parent->parser()->s(), parent->parser()->a3_real(),
          parent->parser()->a3_imag(), parent->parser()->m(), ifactory) {
}

IR3 eigenmode_castor_a::covariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double phi = position[IR3::w];
  double chi = metric_->reduce_chi(position[IR3::v]);
  std::complex<double> factor =
      native_factor_*std::exp(eigenvalue_*time + i_n_tor_*phi);
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
