// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Paulo Rodrigues and Manuel Assunção.

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

// @morphism_helena.cc, this file is part of ::gyronimo::

#include <gyronimo/core/multiroot.hh>
#include <gyronimo/metrics/morphism_helena.hh>

#include <cmath>
#include <numbers>

namespace gyronimo {

morphism_helena::morphism_helena(
    const parser_helena* p, const interpolator2d_factory* ifactory)
    : parser_(p), R_(nullptr), z_(nullptr) {
  double Rgeo = p->rgeo();
  double a = p->eps() * Rgeo;
  dblock_adapter s_range(p->s()), chi_range(p->chi());
  R_ = ifactory->interpolate_data(
      s_range, chi_range,
      dblock_adapter(parser_helena::narray_type(a * p->x() + Rgeo)));
  z_ = ifactory->interpolate_data(
      s_range, chi_range,
      dblock_adapter(parser_helena::narray_type(a * p->y() + 0.0)));
}
morphism_helena::~morphism_helena() {
  if (R_) delete R_;
  if (z_) delete z_;
}
IR3 morphism_helena::operator()(const IR3& q) const {
  double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]), phi = q[IR3::w];
  double R = (*R_)(s, chi);
  return {R * std::cos(phi), -R * std::sin(phi), (*z_)(s, chi)};
}
IR3 morphism_helena::inverse(const IR3& X) const {
  double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
  double R = std::sqrt(x * x + y * y);
  multiroot root_finder(gsl_multiroot_fsolver_hybrids, 1.0e-12, 75);
  using IR2 = std::array<double, 2>;
  IR2 guess = {0.5, std::atan2(z, R - parser_->rmag())};
  std::function<IR2(const IR2&)> zero_function = [&](const IR2& args) {
    auto [s, chi] = reflection_past_axis(args[0], args[1]);
    return IR2({(*R_)(s, chi) - R, (*z_)(s, chi) - z});
  };
  IR2 roots = root_finder(zero_function, guess);
  auto [s, chi] = reflection_past_axis(roots[0], roots[1]);
  return {s, chi, std::atan2(-y, x)};
}
IR3 morphism_helena::translation(const IR3& q, const IR3& delta) const {
  IR3 X = (*this)(q);
  double x = X[IR3::u] + delta[IR3::u];
  double y = X[IR3::v] + delta[IR3::v];
  double z = X[IR3::w] + delta[IR3::w];
  double R = std::sqrt(x * x + y * y);
  multiroot root_finder(gsl_multiroot_fsolver_hybrids, 1.0e-12, 75);
  using IR2 = std::array<double, 2>;
  IR2 guess = {q[IR3::u], q[IR3::v]};
  std::function<IR2(const IR2&)> zero_function = [&](const IR2& args) {
    auto [s, chi] = reflection_past_axis(args[0], args[1]);
    return IR2 {(*R_)(s, chi) - R, (*z_)(s, chi) - z};
  };
  IR2 roots = root_finder(zero_function, guess);
  auto [s, chi] = reflection_past_axis(roots[0], roots[1]);
  return {s, chi, std::atan2(-y, x)};
}
dIR3 morphism_helena::del(const IR3& q) const {
  double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]), phi = q[IR3::w];
  double R = (*R_)(s, chi);
  double Ru = R_->partial_u(s, chi), Rv = R_->partial_v(s, chi);
  double cos = std::cos(phi), sin = std::sin(phi);
  return {Ru * cos, Rv * cos, -R * sin, -Ru * sin, -Rv * sin,
      -R * cos, z_->partial_u(s, chi), z_->partial_v(s, chi), 0.0};
}
ddIR3 morphism_helena::ddel(const IR3& q) const {
  double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]), phi = q[IR3::w];
  double R = (*R_)(s, chi);
  double Ru = R_->partial_u(s, chi), Rv = R_->partial_v(s, chi);
  double Ruu = R_->partial2_uu(s, chi), Ruv = R_->partial2_uv(s, chi),
      Rvv = R_->partial2_vv(s, chi);
  double Zu = z_->partial_u(s, chi), Zv = z_->partial_v(s, chi);
  double Zuu = z_->partial2_uu(s, chi), Zuv = z_->partial2_uv(s, chi),
      Zvv = z_->partial2_vv(s, chi);
  double cos = std::cos(phi), sin = std::sin(phi);
  return {Ruu * cos, Ruv * cos, -Ru * sin, Rvv * cos,
      -Rv * sin, -R * cos, -Ruu * sin, -Ruv * sin, -Ru * cos,
      -Rvv * sin, -Rv * cos, R * sin, Zuu, Zuv, 0, Zvv, 0, 0};
}
double morphism_helena::jacobian(const IR3& q) const {
  double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]);
  double R = (*R_)(s, chi);
  double Ru = R_->partial_u(s, chi), Rv = R_->partial_v(s, chi);
  double Zu = z_->partial_u(s, chi), Zv = z_->partial_v(s, chi);
  return R * (Ru * Zv - Rv * Zu);
}
std::pair<double, double> morphism_helena::reflection_past_axis(
    double s, double chi) const {
  if (s < 0) return {-s, parser_->reduce_chi(chi + std::numbers::pi)};
  else return {s, parser_->reduce_chi(chi)};
}

}  // end namespace gyronimo
