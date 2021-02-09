// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @eigenmode_castor_a.cc

#include <ranges>
#include <gyronimo/fields/eigenmode_castor_a.hh>

namespace gyronimo{

eigenmode_castor_a::eigenmode_castor_a(
    double m_factor, double t_factor,
    const parser_castor *p, const metric_helena *g,
    const interpolator1d_factory* ifactory)
    : IR3field(m_factor, t_factor, g),
      norm_factor_(1.0), parser_(p), metric_(g),
      w_(p->w_real(), p->w_imag()), i_n_tor_(0.0, p->n_tor()),
      tildeA1_(p->s(), p->a1_real(), p->a1_imag(), p->m(), ifactory),
      tildeA2_(p->s(), p->a2_real(), p->a2_imag(), p->m(), ifactory),
      tildeA3_(p->s(), p->a3_real(), p->a3_imag(), p->m(), ifactory) {
  norm_factor_ = 1.0/std::ranges::max(
      parser_->s() | std::views::transform(
          [this](double s){return this->magnitude({s, 0, 0}, 0);}));
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
