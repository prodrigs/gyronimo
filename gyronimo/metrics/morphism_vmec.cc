// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Manuel Assunção.

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

#include <numbers>

namespace gyronimo {

morphism_vmec::morphism_vmec(
    const parser_vmec* p, const interpolator1d_factory* ifactory)
    : parser_(p), xm_(p->xm()), xn_(p->xn()) {
  dblock_adapter s_range(p->radius());
  Rmnc_ = new interpolator1d*[xm_.size()];
  Zmns_ = new interpolator1d*[xm_.size()];
  for (size_t i = 0; i < xm_.size(); ++i) {
    std::slice s_cut(i, s_range.size(), xm_.size());
    narray_type rmnc_i = (p->rmnc())[s_cut];
    Rmnc_[i] = ifactory->interpolate_data(s_range, dblock_adapter(rmnc_i));
    narray_type zmns_i = (p->zmns())[s_cut];
    Zmns_[i] = ifactory->interpolate_data(s_range, dblock_adapter(zmns_i));
  }
}
morphism_vmec::~morphism_vmec() {
  if (Rmnc_) {
    for (size_t i = 0; i < xm_.size(); ++i)
      if (Rmnc_[i]) delete Rmnc_[i];
    delete Rmnc_;
  }
  if (Zmns_) {
    for (size_t i = 0; i < xm_.size(); ++i)
      if (Zmns_[i]) delete Zmns_[i];
    delete Zmns_;
  }
}
IR3 morphism_vmec::operator()(const IR3& q) const {
  double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
  double R = 0.0, Z = 0.0;

#pragma omp parallel for reduction(+ : R, Z)
  for (size_t i = 0; i < xm_.size(); ++i) {
    double m = xm_[i], n = xn_[i];
    double angle_mn = m * theta - n * zeta;
    double cosmn = std::cos(angle_mn), sinmn = std::sin(angle_mn);
    R += (*Rmnc_[i])(s)*cosmn; // assuming stellarator symmetry
    Z += (*Zmns_[i])(s)*sinmn;
  }
  return {R * std::cos(zeta), R * std::sin(zeta), Z};
}
IR3 morphism_vmec::inverse(const IR3& X) const {
  typedef std::array<double, 2> IR2;
  double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
  double r = std::sqrt(x * x + y * y), zeta = std::atan2(y, x);

  std::function<IR2(const IR2&)> zero_function = [&](const IR2& args) {
    auto [s, theta] = reflection_past_axis(args[0], args[1]);
    auto [R, Z] = get_rz(IR3({s, zeta, theta}));
    return IR2({R - r, Z - z});
  };
  auto [R0, Z0] = get_rz({0, zeta, 0});
  IR2 guess = {0.5, std::atan2(z - Z0, r - R0)};
  IR2 roots = multiroot(1.0e-13, 100)(zero_function, guess);
  auto [s, theta] = reflection_past_axis(roots[0], roots[1]);
  return {s, zeta, theta};
}
IR3 morphism_vmec::translation(const IR3& q, const IR3& delta) const {
  typedef std::array<double, 2> IR2;
  IR3 X = (*this)(q);
  double Xt = X[IR3::u] + delta[IR3::u];
  double Yt = X[IR3::v] + delta[IR3::v];
  double Zt = X[IR3::w] + delta[IR3::w];
  double Rt = std::sqrt(Xt * Xt + Yt * Yt);
  double zeta = std::atan2(Yt, Xt);
  std::function<IR2(const IR2&)> zero_function = [&](const IR2& args) {
    auto [s, theta] = reflection_past_axis(args[0], args[1]);
    auto [R, Z] = get_rz(IR3({s, zeta, theta}));
    return IR2({R - Rt, Z - Zt});
  };
  IR2 guess = {q[IR3::u], q[IR3::w]};
  IR2 roots = multiroot(1.0e-12, 100)(zero_function, guess);
  auto [s, theta] = reflection_past_axis(roots[0], roots[1]);
  return {s, zeta, theta};
}
dIR3 morphism_vmec::del(const IR3& q) const {
  double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
  double R = 0.0, dR_ds = 0.0, dR_dtheta = 0.0, dR_dzeta = 0.0;
  double dZ_ds = 0.0, dZ_dtheta = 0.0, dZ_dzeta = 0.0;
  double sn_zeta = std::sin(zeta);
  double cn_zeta = std::cos(zeta);
#pragma omp parallel for reduction(+: R, dR_ds, dR_dtheta, dR_dzeta, dZ_ds, dZ_dtheta, dZ_dzeta)
  for (size_t i = 0; i < xm_.size(); ++i) {
    double m = xm_[i];
    double n = xn_[i];
    double angle_mn = m * theta - n * zeta;
    double cosmn = std::cos(angle_mn), sinmn = std::sin(angle_mn);
    double rmnc_i = (*Rmnc_[i])(s), zmns_i = (*Zmns_[i])(s);
    R += rmnc_i * cosmn; // assuming stellarator symmetry
    dR_ds += (*Rmnc_[i]).derivative(s) * cosmn;
    dR_dtheta += m * rmnc_i * sinmn;
    dR_dzeta += n * rmnc_i * sinmn;
    dZ_ds += (*Zmns_[i]).derivative(s) * sinmn;
    dZ_dtheta += m * zmns_i * cosmn;
    dZ_dzeta += n * zmns_i * cosmn;
  }
  dR_dtheta = -dR_dtheta;
  dZ_dzeta = -dZ_dzeta;
  return {dR_ds * cn_zeta, dR_dzeta * cn_zeta - R * sn_zeta,
      dR_dtheta * cn_zeta, dR_ds * sn_zeta, dR_dzeta * sn_zeta + R * cn_zeta,
      dR_dtheta * sn_zeta, dZ_ds, dZ_dzeta, dZ_dtheta};
}
ddIR3 morphism_vmec::ddel(const IR3& q) const {
  double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
  double cn = std::cos(zeta), sn = std::sin(zeta);
  double R = 0.0, dR_ds = 0.0, dR_dtheta = 0.0, dR_dzeta = 0.0;
  double d2R_ds2 = 0.0, d2R_dsdtheta = 0.0, d2R_dsdzeta = 0.0;
  double d2R_dzeta2 = 0.0, d2R_dzetadtheta = 0.0, d2R_dtheta2 = 0.0;
  double dZ_ds = 0.0, dZ_dtheta = 0.0, dZ_dzeta = 0.0;
  double d2Z_ds2 = 0.0, d2Z_dsdtheta = 0.0, d2Z_dsdzeta = 0.0;
  double d2Z_dzeta2 = 0.0, d2Z_dzetadtheta = 0.0, d2Z_dtheta2 = 0.0;

#pragma omp parallel for reduction(+: R, dR_ds, dR_dtheta, dR_dzeta, d2R_ds2, d2R_dsdtheta, d2R_dsdzeta, d2R_dtheta2, d2R_dthetadzeta, d2R_dzeta2, Z, dZ_ds ,dZ_dtheta, dZ_dzeta, d2Z_ds2, d2Z_dsdtheta, d2Z_dsdzeta, d2Z_dtheta2, d2Z_dthetadzeta, d2Z_dzeta2)
  for (size_t i = 0; i < xm_.size(); ++i) {
    double m = xm_[i], n = xn_[i];
    double cosmn = std::cos(m * theta - n * zeta);
    double sinmn = std::sin(m * theta - n * zeta);
    double rmnc_i = (*Rmnc_[i])(s), zmns_i = (*Zmns_[i])(s);
    double d_rmnc_i = (*Rmnc_[i]).derivative(s);
    double d_zmns_i = (*Zmns_[i]).derivative(s);
    double d2_rmnc_i = (*Rmnc_[i]).derivative2(s);
    double d2_zmns_i = (*Zmns_[i]).derivative2(s);
    R += rmnc_i * cosmn; // assuming stellarator symmetry
    dR_ds += d_rmnc_i * cosmn;
    dR_dzeta += n * rmnc_i * sinmn;
    dR_dtheta -= m * rmnc_i * sinmn;
    d2R_ds2 += d2_rmnc_i * cosmn;
    d2R_dsdzeta += n * d_rmnc_i * sinmn;
    d2R_dsdtheta -= m * d_rmnc_i * sinmn;
    d2R_dzeta2 -= n * n * rmnc_i * cosmn;
    d2R_dzetadtheta += m * n * rmnc_i * cosmn;
    d2R_dtheta2 -= m * m * rmnc_i * cosmn;
    dZ_ds += d_zmns_i * sinmn;
    dZ_dzeta -= n * zmns_i * cosmn;
    dZ_dtheta += m * zmns_i * cosmn;
    d2Z_ds2 += d2_zmns_i * sinmn;
    d2Z_dsdzeta -= n * d_zmns_i * cosmn;
    d2Z_dsdtheta += m * d_zmns_i * cosmn;
    d2Z_dzeta2 -= n * n * zmns_i * sinmn;
    d2Z_dzetadtheta += m * n * zmns_i * sinmn;
    d2Z_dtheta2 -= m * m * zmns_i * sinmn;
  }
  return {d2R_ds2 * cn, d2R_dsdzeta * cn - dR_ds * sn, d2R_dsdtheta * cn,
      (d2R_dzeta2 - R) * cn - 2 * dR_dzeta * sn,
      d2R_dzetadtheta * cn - dR_dtheta * sn, d2R_dtheta2 * cn, d2R_ds2 * sn,
      d2R_dsdzeta * sn + dR_ds * cn, d2R_dsdtheta * sn,
      (d2R_dzeta2 - R) * sn + 2 * dR_dzeta * cn,
      d2R_dzetadtheta * sn + dR_dtheta * cn,
      d2R_dtheta2 * sn, d2Z_ds2, d2Z_dsdzeta,
      d2Z_dsdtheta, d2Z_dzeta2, d2Z_dzetadtheta, d2Z_dtheta2};
}
double morphism_vmec::jacobian(const IR3& q) const {
  double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
  double R = 0.0, dR_ds = 0.0, dR_dtheta = 0.0;
  double dZ_ds = 0.0, dZ_dtheta = 0.0;
  double sn_zeta = std::sin(zeta), cn_zeta = std::cos(zeta);

#pragma omp parallel for reduction(+: R, Z, dR_ds, dR_dtheta, dR_dzeta, dZ_ds, dZ_dtheta, dZ_dzeta)
  for (size_t i = 0; i < xm_.size(); ++i) {
    double m = xm_[i], n = xn_[i];
    double angle_mn = m * theta - n * zeta;
    double cosmn = std::cos(angle_mn), sinmn = std::sin(angle_mn);
    double rmnc_i = (*Rmnc_[i])(s), zmns_i = (*Zmns_[i])(s);
    R += rmnc_i * cosmn; // assuming stellarator symmetry
    dR_ds += (*Rmnc_[i]).derivative(s) * cosmn;
    dR_dtheta += m * rmnc_i * sinmn;
    dZ_ds += (*Zmns_[i]).derivative(s) * sinmn;
    dZ_dtheta += m * zmns_i * cosmn;
  }
  dR_dtheta = -dR_dtheta;
  return R * (dR_ds * dZ_dtheta - dR_dtheta * dZ_ds);
}
std::pair<double, double> morphism_vmec::get_rz(const IR3& position) const {
  double u = position[gyronimo::IR3::u];
  double v = position[gyronimo::IR3::v];
  double w = position[gyronimo::IR3::w];
  double R = 0.0, Z = 0.0;

#pragma omp parallel for reduction(+ : R, Z)
  for (size_t i = 0; i < xm_.size(); ++i) {
    double m = xm_[i], n = xn_[i];
    R += (*Rmnc_[i])(u)*std::cos(m * w - n * v);
    Z += (*Zmns_[i])(u)*std::sin(m * w - n * v);
  }
  return {R, Z};
}
std::pair<double, double> morphism_vmec::reflection_past_axis(
    double s, double theta) const {
  if (s < 0)
    return {-s, theta + std::numbers::pi};
  else return {s, theta};
}

}  // end namespace gyronimo
