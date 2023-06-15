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

// @equilibrium_vmec.cc, this file is part of ::gyronimo::

#include <cmath>
#include <gyronimo/core/dblock.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>

namespace gyronimo {

equilibrium_vmec::equilibrium_vmec(
    const metric_vmec* g, const interpolator1d_factory* ifactory)
    : IR3field_c1(std::abs(g->parser()->B0()), 1.0, g),
      harmonics_(g->parser()->mnmax_nyq()), metric_(g),
      m_(g->parser()->xm_nyq()), n_(g->parser()->xn_nyq()),
      btheta_mn_(g->parser()->mnmax_nyq()),
      bzeta_mn_(g->parser()->mnmax_nyq()) {
  const parser_vmec* p = metric_->parser();
  dblock_adapter sgrid(p->sgrid());
#pragma omp parallel for
  for (size_t i = 0; i < harmonics_; i++) {
    std::slice mask(i, sgrid.size(), harmonics_);
    narray_type bzeta_mn_data = (p->bsupvmnc())[mask] / this->m_factor();
    narray_type btheta_mn_data = (p->bsupumnc())[mask] / this->m_factor();
    bzeta_mn_[i] = std::move(std::unique_ptr<interpolator1d>(
        ifactory->interpolate_data(sgrid, dblock_adapter(bzeta_mn_data))));
    btheta_mn_[i] = std::move(std::unique_ptr<interpolator1d>(
        ifactory->interpolate_data(sgrid, dblock_adapter(btheta_mn_data))));
  };
}

IR3 equilibrium_vmec::contravariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  double b_theta = 0.0, b_zeta = 0.0;
  auto cis_mn = equilibrium_vmec::cached_cis(theta, zeta);
#pragma omp parallel for reduction(+ : b_zeta, b_theta)
  for (size_t i = 0; i < harmonics_; i++) {
    double cos_mn_i = std::real(cis_mn[i]);
    b_zeta += (*bzeta_mn_[i])(s) * cos_mn_i;
    b_theta += (*btheta_mn_[i])(s) * cos_mn_i;
  };
  return {0.0, b_zeta, b_theta};
}

dIR3 equilibrium_vmec::del_contravariant(
    const IR3& position, double time) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  auto cis_mn = equilibrium_vmec::cached_cis(theta, zeta);
  double db_theta_ds = 0.0, db_theta_dtheta = 0.0, db_theta_dzeta = 0.0;
  double db_zeta_ds = 0.0, db_zeta_dtheta = 0.0, db_zeta_dzeta = 0.0;
#pragma omp parallel for reduction(+: db_theta_ds, db_theta_dtheta, db_theta_dzeta, db_zeta_ds, db_zeta_dtheta, db_zeta_dzeta)
  for (size_t i = 0; i < harmonics_; i++) {
    double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
    double bzeta_mn_i =(*bzeta_mn_[i])(s), btheta_mn_i = (*btheta_mn_[i])(s);
    db_theta_ds += (*btheta_mn_[i]).derivative(s) * cos_mn_i;
    db_theta_dtheta -= m_[i] * btheta_mn_i * sin_mn_i;
    db_theta_dzeta += n_[i] * btheta_mn_i * sin_mn_i;
    db_zeta_ds += (*bzeta_mn_[i]).derivative(s) * cos_mn_i;
    db_zeta_dtheta -= m_[i] * bzeta_mn_i * sin_mn_i;
    db_zeta_dzeta += n_[i] * bzeta_mn_i * sin_mn_i;
  };
  return {
      0.0,
      0.0,
      0.0,
      db_zeta_ds,
      db_zeta_dzeta,
      db_zeta_dtheta,
      db_theta_ds,
      db_theta_dzeta,
      db_theta_dtheta};
}

//! Cached evaluation of trigonometric coefficients for internal Fourier series.
const equilibrium_vmec::cis_container_t& equilibrium_vmec::cached_cis(
    double theta, double zeta) const {
  thread_local double cached_theta = -1e6, cached_zeta = -1e6;
  thread_local cis_container_t cis_mn(harmonics_);
  if (theta != cached_theta || zeta != cached_zeta) {
#pragma omp parallel for
    for(size_t i = 0;i < harmonics_;i++) {
      double angle_mn = m_[i] * theta - n_[i] * zeta;
      cis_mn[i] = {std::cos(angle_mn), std::sin(angle_mn)};
    }
    cached_theta = theta;
    cached_zeta = zeta;
  }
  return cis_mn;
}

}  // end namespace gyronimo.
