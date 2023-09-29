// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Jorge Ferreira and Paulo Rodrigues.

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

// @equilibrium_vmec.cc, this file is part of ::gyronimo::

#include <gyronimo/core/dblock.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>

#include <cmath>
#include <numeric>

namespace gyronimo {

equilibrium_vmec::equilibrium_vmec(
    const metric_vmec* g, const interpolator1d_factory* ifactory)
    : IR3field_c1(std::abs(g->my_parser()->B0()), 1.0, g),
      metric_(g), parser_(g->my_parser()), harmonics_(parser_->mnmax_nyq()),
      m_(parser_->xm_nyq()), n_(parser_->xn_nyq()), index_(harmonics_),
      btheta_mn_(parser_->mnmax_nyq()), bzeta_mn_(parser_->mnmax_nyq()) {
  std::iota(index_.begin(), index_.end(), 0);
  this->build_interpolator_array(
      bzeta_mn_, parser_->bsupvmnc() / this->m_factor(), ifactory);
  this->build_interpolator_array(
      btheta_mn_, parser_->bsupumnc() / this->m_factor(), ifactory);
}

IR3 equilibrium_vmec::contravariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  const auto& cis_mn = equilibrium_vmec::cached_cis(theta, zeta);
  auto out = std::transform_reduce(
      index_.begin(), index_.end(), auxiliar1_t {0, 0}, std::plus<>(),
      [&](size_t i) -> auxiliar1_t {
        double cos_mn = std::real(cis_mn[i]);
        return {(*bzeta_mn_[i])(s)*cos_mn, (*btheta_mn_[i])(s)*cos_mn};
      });
  return {0, out.bzeta, out.btheta};
}

dIR3 equilibrium_vmec::del_contravariant(
    const IR3& position, double time) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  const auto& cis_mn = equilibrium_vmec::cached_cis(theta, zeta);
  auto out = std::transform_reduce(
      index_.begin(), index_.end(), auxiliar2_t {0, 0, 0, 0, 0, 0},
      std::plus<>(), [&](size_t i) -> auxiliar2_t {
        double bzeta_mn_i = (*bzeta_mn_[i])(s);
        double btheta_mn_i = (*btheta_mn_[i])(s);
        double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
        return {
            (*bzeta_mn_[i]).derivative(s) * cos_mn_i,
            n_[i] * bzeta_mn_i * sin_mn_i,
            -m_[i] * bzeta_mn_i * sin_mn_i,
            (*btheta_mn_[i]).derivative(s) * cos_mn_i,
            n_[i] * btheta_mn_i * sin_mn_i,
            -m_[i] * btheta_mn_i * sin_mn_i};
      });
  return {0.0, 0.0, 0.0,
      out.dbzetadu, out.dbzetadv, out.dbzetadw,
      out.dbthetadu, out.dbthetadv, out.dbthetadw};
}

void equilibrium_vmec::build_interpolator_array(
    std::vector<std::unique_ptr<interpolator1d>>& interpolator_array,
    const narray_type& samples_array, const interpolator1d_factory* ifactory) {
  dblock_adapter sgrid(parser_->sgrid());
  std::transform(
      index_.begin(), index_.end(), interpolator_array.begin(), [&](size_t i) {
        std::slice mask_i(i, sgrid.size(), harmonics_);
        narray_type data = samples_array[mask_i];
        return std::move(std::unique_ptr<interpolator1d>(
            ifactory->interpolate_data(sgrid, dblock_adapter(data))));
      });
}

//! Cached evaluation of trigonometric coefficients for internal Fourier series.
const equilibrium_vmec::cis_container_t& equilibrium_vmec::cached_cis(
    double theta, double zeta) const {
  thread_local double cached_theta = -1e6, cached_zeta = -1e6;
  thread_local cis_container_t cis_mn(harmonics_);
  if (theta != cached_theta || zeta != cached_zeta) {
    std::transform(
        index_.begin(), index_.end(), cis_mn.begin(),
        [&](size_t i) -> cis_container_t::value_type {
          double angle_mn = m_[i] * theta - n_[i] * zeta;
          return {std::cos(angle_mn), std::sin(angle_mn)};
        });
    cached_theta = theta;
    cached_zeta = zeta;
  }
  return cis_mn;
}

}  // end namespace gyronimo.
