// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues and Manuel Assunção.

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

// @metric_helena.cc, this file is part of ::gyronimo::

#include <gyronimo/core/transpose.hh>
#include <gyronimo/metrics/metric_helena.hh>

#include <numbers>

namespace gyronimo {

metric_helena::metric_helena(
    const morphism_helena* morph, const interpolator2d_factory* ifactory)
    : metric_connected(morph), parser_(morph->parser()), R0_(parser_->rmag()),
      squaredR0_(parser_->rmag() * parser_->rmag()), guu_(nullptr),
      guv_(nullptr), gvv_(nullptr), gww_(nullptr) {
  dblock_adapter s_range(parser_->s()), chi_range(parser_->chi());
  guu_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(parser_->covariant_g11()));
  guv_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(parser_->covariant_g12()));
  gvv_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(parser_->covariant_g22()));
  gww_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(parser_->covariant_g33()));
}
metric_helena::~metric_helena() {
  if (guu_) delete guu_;
  if (guv_) delete guv_;
  if (gvv_) delete gvv_;
  if (gww_) delete gww_;
}
SM3 metric_helena::operator()(const IR3& q) const {
  double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]);
  return {
      squaredR0_ * (*guu_)(s, chi), squaredR0_ * (*guv_)(s, chi), 0,
      squaredR0_ * (*gvv_)(s, chi), 0, squaredR0_ * (*gww_)(s, chi)};
}
dSM3 metric_helena::del(const IR3& q) const {
  double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]);
  return {
      squaredR0_ * (*guu_).partial_u(s, chi),
      squaredR0_ * (*guu_).partial_v(s, chi), 0,
      squaredR0_ * (*guv_).partial_u(s, chi),
      squaredR0_ * (*guv_).partial_v(s, chi), 0, 0, 0, 0,
      squaredR0_ * (*gvv_).partial_u(s, chi),
      squaredR0_ * (*gvv_).partial_v(s, chi), 0, 0, 0, 0,
      squaredR0_ * (*gww_).partial_u(s, chi),
      squaredR0_ * (*gww_).partial_v(s, chi), 0};
}

}  // end namespace gyronimo
