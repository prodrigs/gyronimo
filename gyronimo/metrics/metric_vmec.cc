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

namespace gyronimo {

metric_vmec::metric_vmec(
    const parser_vmec* p, const interpolator1d_factory* ifactory)
    : parser_(p), harmonics_(p->mnmax()), m_(p->xm()), n_(p->xn()),
    r_mn_(p->mnmax()), z_mn_(p->mnmax()) {
  dblock_adapter sgrid(p->sgrid());
#pragma omp parallel for
  for (size_t i = 0; i < harmonics_; i++) {
    std::slice mask(i, sgrid.size(), harmonics_);
    narray_type r_mn_data = (p->rmnc())[mask], z_mn_data = (p->zmns())[mask];
    r_mn_[i] = std::move(std::unique_ptr<interpolator1d>(
        ifactory->interpolate_data(sgrid, dblock_adapter(r_mn_data))));
    z_mn_[i] = std::move(std::unique_ptr<interpolator1d>(
        ifactory->interpolate_data(sgrid, dblock_adapter(z_mn_data))));
  }
}

SM3 metric_vmec::operator()(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  auto cis_mn = metric_vmec::cached_cis(theta, zeta);
  double r = 0.0, dr_ds = 0.0, dr_dtheta = 0.0, dr_dzeta = 0.0;
  double dz_ds = 0.0, dz_dtheta = 0.0, dz_dzeta = 0.0;
#pragma omp parallel for reduction(+: r, z, dr_ds, dr_dtheta, dr_dzeta, dz_ds, dz_dtheta, dz_dzeta)
  for (size_t i = 0; i < harmonics_; i++) {
    double r_mn_i = (*r_mn_[i])(s), z_mn_i = (*z_mn_[i])(s);
    double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
    r += r_mn_i * cos_mn_i;
    dr_ds += (*r_mn_[i]).derivative(s) * cos_mn_i;
    dr_dtheta -= m_[i] * r_mn_i * sin_mn_i;
    dr_dzeta += n_[i] * r_mn_i * sin_mn_i;
    dz_ds += (*z_mn_[i]).derivative(s) * sin_mn_i;
    dz_dtheta += m_[i] * z_mn_i * cos_mn_i;
    dz_dzeta -= n_[i] * z_mn_i * cos_mn_i;
  };
  return {
      dr_ds * dr_ds + dz_ds * dz_ds,  // g_uu
      dr_ds * dr_dzeta + dz_ds * dz_dzeta,  // g_uw
      dr_ds * dr_dtheta + dz_ds * dz_dtheta,  // g_uv
      r * r + dr_dzeta * dr_dzeta + dz_dzeta * dz_dzeta,  // g_vv
      dr_dtheta * dr_dzeta + dz_dtheta * dz_dzeta,  // g_vw
      dr_dtheta * dr_dtheta + dz_dtheta * dz_dtheta  // g_ww
  };
}

dSM3 metric_vmec::del(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  auto cis_mn = metric_vmec::cached_cis(theta, zeta);
  double r = 0.0, z = 0.0;
  double dr_ds = 0.0, dr_dtheta = 0.0, dr_dzeta = 0.0;
  double d2r_ds2 = 0.0, d2r_dsdtheta = 0.0, d2r_dsdzeta = 0.0;
  double d2r_dthetads = 0.0, d2r_dtheta2 = 0.0, d2r_dthetadzeta = 0.0;
  double d2r_dzetads = 0.0, d2r_dzetadtheta = 0.0, d2r_dzeta2 = 0.0;
  double dz_ds = 0.0, dz_dtheta = 0.0, dz_dzeta = 0.0;
  double d2z_ds2 = 0.0, d2z_dsdtheta = 0.0, d2z_dsdzeta = 0.0;
  double d2z_dthetads = 0.0, d2z_dtheta2 = 0.0, d2z_dthetadzeta = 0.0;
  double d2z_dzetads = 0.0, d2z_dzetadtheta = 0.0, d2z_dzeta2 = 0.0;
#pragma omp parallel for reduction(+: r, dr_ds, dr_dtheta, dr_dzeta, d2r_ds2, d2r_dsdtheta, d2r_dsdzeta, d2r_dtheta2, d2r_dthetadzeta, d2r_dzeta2, z, dz_ds ,dz_dtheta, dz_dzeta, d2z_ds2, d2z_dsdtheta, d2z_dsdzeta, d2z_dtheta2, d2z_dthetadzeta, d2z_dzeta2)
  for (size_t i = 0; i < harmonics_; i++) {
    double r_mn_i = (*r_mn_[i])(s), z_mn_i = (*z_mn_[i])(s);
    double cos_mn_i = std::real(cis_mn[i]), sin_mn_i = std::imag(cis_mn[i]);
    double d_r_mn_i = (*r_mn_[i]).derivative(s);
    double d_z_mn_i = (*z_mn_[i]).derivative(s);
    double d2_r_mn_i = (*r_mn_[i]).derivative2(s);
    double d2_z_mn_i = (*z_mn_[i]).derivative2(s);
    r += r_mn_i * cos_mn_i;
    z += z_mn_i * sin_mn_i;
    dr_ds += d_r_mn_i * cos_mn_i;
    dr_dtheta -= m_[i] * r_mn_i * sin_mn_i;
    dr_dzeta += n_[i] * r_mn_i * sin_mn_i;
    d2r_ds2 += d2_r_mn_i * cos_mn_i;
    d2r_dsdtheta -= m_[i] * d_r_mn_i * sin_mn_i;
    d2r_dsdzeta += n_[i] * d_r_mn_i * sin_mn_i;
    d2r_dtheta2 -= m_[i] * m_[i] * r_mn_i * cos_mn_i;
    d2r_dthetadzeta += m_[i] * n_[i] * r_mn_i * cos_mn_i;
    d2r_dzeta2 -= n_[i] * n_[i] * r_mn_i * cos_mn_i;
    dz_ds += d_z_mn_i * sin_mn_i;
    dz_dtheta += m_[i] * z_mn_i * cos_mn_i;
    dz_dzeta -= n_[i] * z_mn_i * cos_mn_i;
    d2z_ds2 += d2_z_mn_i * sin_mn_i;
    d2z_dsdtheta += m_[i] * d_z_mn_i * cos_mn_i;
    d2z_dsdzeta -= n_[i] * d_z_mn_i * cos_mn_i;
    d2z_dtheta2 -= m_[i] * m_[i] * z_mn_i * sin_mn_i;
    d2z_dthetadzeta += m_[i] * n_[i] * z_mn_i * sin_mn_i;
    d2z_dzeta2 -= n_[i] * n_[i] * z_mn_i * sin_mn_i;
  }
  return {
      2 * (dr_ds * d2r_ds2 + dz_ds * d2z_ds2),
      2 * (dr_ds * d2r_dsdzeta + dz_ds * d2z_dsdzeta),  // d_i g_uu
      2 * (dr_ds * d2r_dsdtheta + dz_ds * d2z_dsdtheta),
      dr_ds * d2r_dsdzeta + dr_dzeta * d2r_ds2 + dz_ds * d2z_dsdzeta +
          dz_dzeta * d2z_ds2,
      dr_ds * d2r_dzeta2 + dr_dzeta * d2r_dsdzeta + dz_ds * d2z_dzeta2 +
          dz_dzeta * d2z_dsdzeta,  // d_i g_uv
      dr_ds * d2r_dthetadzeta + dr_dzeta * d2r_dsdtheta +
          dz_ds * d2z_dthetadzeta + dz_dzeta * d2z_dsdtheta,
      dr_ds * d2r_dsdtheta + dr_dtheta * d2r_ds2 + dz_ds * d2z_dsdtheta +
          dz_dtheta * d2z_ds2,
      dr_ds * d2r_dthetadzeta + dr_dtheta * d2r_dsdzeta +
          dz_ds * d2z_dthetadzeta + dz_dtheta * d2z_dsdzeta,  // d_i g_uw
      dr_ds * d2r_dtheta2 + dr_dtheta * d2r_dsdtheta + dz_ds * d2z_dtheta2 +
          dz_dtheta * d2z_dsdtheta,
      2 * (r * dr_ds + dr_dzeta * d2r_dsdzeta + dz_dzeta * d2z_dsdzeta),
      2 *
          (r * dr_dzeta + dr_dzeta * d2r_dzeta2 +
           dz_dzeta * d2z_dzeta2),  // d_i g_vv
      2 *
          (r * dr_dtheta + dr_dzeta * d2r_dthetadzeta +
           dz_dzeta * d2z_dthetadzeta),
      dr_dtheta * d2r_dsdzeta + dr_dzeta * d2r_dsdtheta +
          dz_dtheta * d2z_dsdzeta + dz_dzeta * d2z_dsdtheta,
      dr_dtheta * d2r_dzeta2 + dr_dzeta * d2r_dthetadzeta +
          dz_dtheta * d2z_dzeta2 + dz_dzeta * d2z_dthetadzeta,  // d_i g_vw
      dr_dtheta * d2r_dthetadzeta + dr_dzeta * d2r_dtheta2 +
          dz_dtheta * d2z_dthetadzeta + dz_dzeta * d2z_dtheta2,
      2 * (dr_dtheta * d2r_dsdtheta + dz_dtheta * d2z_dsdtheta),
      2 * (dr_dtheta * d2r_dthetadzeta + dz_dtheta * d2z_dthetadzeta),  // d_i
                                                                        // g_ww
      2 * (dr_dtheta * d2r_dtheta2 + dz_dtheta * d2z_dtheta2),
  };
}

IR3 metric_vmec::transform2cylindrical(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  auto cis_mn = metric_vmec::cached_cis(theta, zeta);
  double r = 0.0, z = 0.0;
#pragma omp parallel for reduction(+ : r, z)
  for (size_t i = 0; i < harmonics_; i++) {
    double r_mn_i = (*r_mn_[i])(s), z_mn_i = (*z_mn_[i])(s);
    r += r_mn_i * std::real(cis_mn[i]);
    z += z_mn_i * std::imag(cis_mn[i]);
  }
  return {r, zeta, z};
}

//! Cached evaluation of trigonometric coefficients for internal Fourier series.
const metric_vmec::cis_container_t& metric_vmec::cached_cis(
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

}  // end namespace gyronimo
