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

// @mrophism_vmec.cc, this file is part of ::gyronimo::

#include <gyronimo/core/multiroot.hh>
#include <gyronimo/metrics/morphism_vmec.hh>

#include <numeric>

namespace gyronimo {

//! Cached evaluation of trigonometric coefficients for internal Fourier series.
const morphism_vmec::cis_container_t& morphism_vmec::cached_cis(
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

morphism_vmec::morphism_vmec(
    const parser_vmec* p, const interpolator1d_factory* ifactory)
    : parser_(p), harmonics_(p->mnmax()), m_(p->xm()), n_(p->xn()),
      index_(harmonics_), r_mn_(p->mnmax()), z_mn_(p->mnmax()) {
  std::iota(index_.begin(), index_.end(), 0);
  this->build_interpolator_array(r_mn_, parser_->rmnc(), ifactory);
  this->build_interpolator_array(z_mn_, parser_->zmns(), ifactory);
}

void morphism_vmec::build_interpolator_array(
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

IR3 morphism_vmec::inverse(const IR3& X) const {
  double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
  double r = std::sqrt(x * x + y * y), zeta = std::atan2(y, x);
  auto [r_axis, z_axis] = get_rz({0, zeta, 0});
  return this->inverse(X, {0.5, std::atan2(z - z_axis, r - r_axis)});
}

std::pair<double, double> morphism_vmec::get_rz(const IR3& q) const {
  double flux = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
  auto cis_mn = morphism_vmec::cached_cis(theta, zeta);
  auto [r, z] = std::transform_reduce(
      index_.begin(), index_.end(), aux_rz_t {0, 0}, std::plus<>(),
      [&](size_t i) -> aux_rz_t {
        double r_mn_i = (*r_mn_[i])(flux), z_mn_i = (*z_mn_[i])(flux);
        double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
        return {r_mn_i * cos_mn_i, z_mn_i * sin_mn_i};
      });
  return {r, z};
}

IR3 morphism_vmec::inverse(
    const IR3& X, const std::pair<double, double>& guess) const {
  double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
  double r = std::sqrt(x * x + y * y), zeta = std::atan2(y, x);
  using IR2 = std::array<double, 2>;
  std::function<IR2(const IR2&)> zero_function = [&](const IR2& args) -> IR2 {
    auto [flux, theta] = reflection_past_axis(args[0], args[1]);
    auto [r_trial, z_trial] = get_rz({flux, zeta, theta});
    return {r_trial - r, z_trial - z};
  };
  auto roots =
      multiroot(1.0e-12, 100)(zero_function, IR2 {guess.first, guess.second});
  auto [flux, theta] = reflection_past_axis(roots[0], roots[1]);
  return {flux, zeta, theta};
}

dIR3 morphism_vmec::del(const IR3& q) const {
  double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
  auto cis_mn = morphism_vmec::cached_cis(theta, zeta);
  auto a = std::transform_reduce(
      index_.begin(), index_.end(), aux_del_t {0, 0, 0, 0, 0, 0, 0},
      std::plus<>(), [&](size_t i) -> aux_del_t {
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
  double sin_zeta = std::sin(zeta), cos_zeta = std::cos(zeta);
  return {
      a.drdu * cos_zeta, a.drdv * cos_zeta - a.r * sin_zeta, a.drdw * cos_zeta,
      a.drdu * sin_zeta, a.drdv * sin_zeta + a.r * cos_zeta, a.drdw * sin_zeta,
      a.dzdu, a.dzdv, a.dzdw};
}

ddIR3 morphism_vmec::ddel(const IR3& q) const {
  double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
  auto cis_mn = morphism_vmec::cached_cis(theta, zeta);
  auto a = std::transform_reduce(
      index_.begin(), index_.end(),
      aux_ddel_t {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      std::plus<>(), [&](size_t i) -> aux_ddel_t {
        double r_mn_i = (*r_mn_[i])(s), z_mn_i = (*z_mn_[i])(s);
        double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
        double drdu_mn_i = (*r_mn_[i]).derivative(s);
        double dzdu_mn_i = (*z_mn_[i]).derivative(s);
        double d2rdudu_mn_i = (*r_mn_[i]).derivative2(s);
        double d2zdudu_mn_i = (*z_mn_[i]).derivative2(s);
        return {
            r_mn_i * cos_mn_i,  // r_mn_i
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
  double sin_zeta = std::sin(zeta), cos_zeta = std::cos(zeta);
  return {
      a.d2rdudu * cos_zeta, a.d2rdudv * cos_zeta - a.drdu * sin_zeta,
      a.d2rdudw * cos_zeta,
      (a.d2rdvdv - a.r) * cos_zeta - 2 * a.drdv * sin_zeta,
      a.d2rdvdw * cos_zeta - a.drdw * sin_zeta, a.d2rdwdw * cos_zeta,
      a.d2rdudu * sin_zeta, a.d2rdudv * sin_zeta + a.drdu * cos_zeta,
      a.d2rdudw * sin_zeta,
      (a.d2rdvdv - a.r) * sin_zeta + 2 * a.drdv * cos_zeta,
      a.d2rdvdw * sin_zeta + a.drdw * cos_zeta, a.d2rdwdw * sin_zeta,
      a.d2zdudu, a.d2zdudv, a.d2zdudw, a.d2zdvdv, a.d2zdvdw, a.d2zdwdw};
}

}  // end namespace gyronimo
