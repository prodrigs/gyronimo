// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022, 2023 Jorge Ferreira and Paulo Rodrigues.

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

// @metric_vmec.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_vmec.hh>

#include <cmath>
#include <numeric>

namespace gyronimo {

metric_vmec::metric_vmec(
    const parser_vmec* p, const interpolator1d_factory* ifactory)
    : parser_(p), harmonics_(p->mnmax()), m_(p->xm()), n_(p->xn()),
      index_(harmonics_), r_mn_(p->mnmax()), z_mn_(p->mnmax()) {
  std::iota(index_.begin(), index_.end(), 0);
  this->build_interpolator_array(r_mn_, parser_->rmnc(), ifactory);
  this->build_interpolator_array(z_mn_, parser_->zmns(), ifactory);
}

SM3 metric_vmec::operator()(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  auto cis_mn = metric_vmec::cached_cis(theta, zeta);
  auto out = std::transform_reduce(
      index_.begin(), index_.end(), auxiliar1_t {0, 0, 0, 0, 0, 0, 0},
      std::plus<>(), [&](size_t i) -> auxiliar1_t {
        double r_mn_i = (*r_mn_[i])(s), z_mn_i = (*z_mn_[i])(s);
        double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
        return {
            r_mn_i * cos_mn_i,  // r_mn_i
            (*r_mn_[i]).derivative(s) * cos_mn_i,  // drdu_mn_i
            n_[i] * r_mn_i * sin_mn_i,  // drdv_mn_i
            -m_[i] * r_mn_i * sin_mn_i,  // drdw_mn_i
            (*z_mn_[i]).derivative(s) * sin_mn_i,  // dzdu_mn_i
            -n_[i] * z_mn_i * cos_mn_i,  // dzdv_mn_i
            m_[i] * z_mn_i * cos_mn_i  // dzdw_mn_i
        };
      });
  return {
      out.drdu * out.drdu + out.dzdu * out.dzdu,  // g_uu
      out.drdu * out.drdv + out.dzdu * out.dzdv,  // g_uw
      out.drdu * out.drdw + out.dzdu * out.dzdw,  // g_uv
      out.r * out.r + out.drdv * out.drdv + out.dzdv * out.dzdv,  // g_vv
      out.drdw * out.drdv + out.dzdw * out.dzdv,  // g_vw
      out.drdw * out.drdw + out.dzdw * out.dzdw  // g_ww
  };
}

dSM3 metric_vmec::del(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  auto cis_mn = metric_vmec::cached_cis(theta, zeta);
  auto out = std::transform_reduce(
      index_.begin(), index_.end(),
      auxiliar2_t {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      std::plus<>(), [&](size_t i) -> auxiliar2_t {
        double r_mn_i = (*r_mn_[i])(s), z_mn_i = (*z_mn_[i])(s);
        double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
        double drdu_mn_i = (*r_mn_[i]).derivative(s);
        double dzdu_mn_i = (*z_mn_[i]).derivative(s);
        double d2rdudu_mn_i = (*r_mn_[i]).derivative2(s);
        double d2zdudu_mn_i = (*z_mn_[i]).derivative2(s);
        return {
            r_mn_i * cos_mn_i,  // r_mn_i
            z_mn_i * sin_mn_i,  // z_mn_i
            drdu_mn_i * cos_mn_i,  // drdu_mn_i
            n_[i] * r_mn_i * sin_mn_i,  // drdv_mn_i
            -m_[i] * r_mn_i * sin_mn_i,  // drdw_mn_i
            dzdu_mn_i * sin_mn_i,  // dzdu_mn_i
            -n_[i] * z_mn_i * cos_mn_i,  // dzdv_mn_i
            m_[i] * z_mn_i * cos_mn_i,  // dzdw_mn_i
            d2rdudu_mn_i * cos_mn_i,  // d2rdudu_mn_i
            n_[i] * drdu_mn_i * sin_mn_i,  // d2rdudv_mn_i
            -m_[i] * drdu_mn_i * sin_mn_i,  // d2rdudw_mn_i
            -n_[i] * n_[i] * r_mn_i * cos_mn_i,  // d2rdvdv_mn_i
            m_[i] * n_[i] * r_mn_i * cos_mn_i,  // d2rdvdw_mn_i
            -m_[i] * m_[i] * r_mn_i * cos_mn_i,  // d2rdwdw_mn_i
            d2zdudu_mn_i * sin_mn_i,  // dzdudu_mn_i
            -n_[i] * dzdu_mn_i * cos_mn_i,  // dzdudv_mn_i
            m_[i] * dzdu_mn_i * cos_mn_i,  // dzdudw_mn_i
            -n_[i] * n_[i] * z_mn_i * sin_mn_i,  // dzdvdv_mn_i
            m_[i] * n_[i] * z_mn_i * sin_mn_i,  // dzdvdw_mn_i
            -m_[i] * m_[i] * z_mn_i * sin_mn_i  // dzdwdw_mn_i
        };
      });
  return {
      2 * (out.drdu * out.d2rdudu + out.dzdu * out.d2zdudu),
      2 * (out.drdu * out.d2rdudv + out.dzdu * out.d2zdudv),
      2 * (out.drdu * out.d2rdudw + out.dzdu * out.d2zdudw),
      out.drdu * out.d2rdudv + out.drdv * out.d2rdudu + out.dzdu * out.d2zdudv +
          out.dzdv * out.d2zdudu,
      out.drdu * out.d2rdvdv + out.drdv * out.d2rdudv + out.dzdu * out.d2zdvdv +
          out.dzdv * out.d2zdudv,
      out.drdu * out.d2rdvdw + out.drdv * out.d2rdudw + out.dzdu * out.d2zdvdw +
          out.dzdv * out.d2zdudw,
      out.drdu * out.d2rdudw + out.drdw * out.d2rdudu + out.dzdu * out.d2zdudw +
          out.dzdw * out.d2zdudu,
      out.drdu * out.d2rdvdw + out.drdw * out.d2rdudv + out.dzdu * out.d2zdvdw +
          out.dzdw * out.d2zdudv,
      out.drdu * out.d2rdwdw + out.drdw * out.d2rdudw + out.dzdu * out.d2zdwdw +
          out.dzdw * out.d2zdudw,
      2 * (out.r * out.drdu + out.drdv * out.d2rdudv + out.dzdv * out.d2zdudv),
      2 * (out.r * out.drdv + out.drdv * out.d2rdvdv + out.dzdv * out.d2zdvdv),
      2 * (out.r * out.drdw + out.drdv * out.d2rdvdw + out.dzdv * out.d2zdvdw),
      out.drdw * out.d2rdudv + out.drdv * out.d2rdudw + out.dzdw * out.d2zdudv +
          out.dzdv * out.d2zdudw,
      out.drdw * out.d2rdvdv + out.drdv * out.d2rdvdw + out.dzdw * out.d2zdvdv +
          out.dzdv * out.d2zdvdw,
      out.drdw * out.d2rdvdw + out.drdv * out.d2rdwdw + out.dzdw * out.d2zdvdw +
          out.dzdv * out.d2zdwdw,
      2 * (out.drdw * out.d2rdudw + out.dzdw * out.d2zdudw),
      2 * (out.drdw * out.d2rdvdw + out.dzdw * out.d2zdvdw),
      2 * (out.drdw * out.d2rdwdw + out.dzdw * out.d2zdwdw),
  };
}

IR3 metric_vmec::transform2cylindrical(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  auto cis_mn = metric_vmec::cached_cis(theta, zeta);
  double r = 0.0, z = 0.0;
  for (size_t i = 0; i < harmonics_; i++) {
    double r_mn_i = (*r_mn_[i])(s), z_mn_i = (*z_mn_[i])(s);
    r += r_mn_i * std::real(cis_mn[i]);
    z += z_mn_i * std::imag(cis_mn[i]);
  }
  return {r, zeta, z};
}

void metric_vmec::build_interpolator_array(
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
const metric_vmec::cis_container_t& metric_vmec::cached_cis(
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

}  // end namespace gyronimo
