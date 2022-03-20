// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2022 Rogerio Jorge and Paulo Rodrigues.

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

// @equilibrium_stellnaqs.cc, this file is part of ::gyronimo::

#include <cmath>
#include <numbers>
#include <gyronimo/fields/equilibrium_stellnaqs.hh>

namespace gyronimo{

equilibrium_stellnaqs::equilibrium_stellnaqs(
    const metric_stellnaqs *g,
    double axis_field, double axis_length, double axis_iota)
    : IR3field_c1(axis_field, 1.0, g), metric_(g), axis_length_(axis_length) {
  length_factor_ = 2.0*std::numbers::pi/axis_length;
  iota_factor_ = (axis_iota - g->field_periods())*length_factor_;
}
IR3 equilibrium_stellnaqs::contravariant(
    const IR3& position, double time) const {
  double r = position[IR3::u], theta = position[IR3::v];
  double phi = metric_->reduce_phi(position[IR3::w]);
  double k = (*metric_->curvature())(phi);
  double Bw = length_factor_*(1.0 + 2.0*k*r*std::cos(theta));
  double Bv = iota_factor_*Bw;
  return {0.0, Bv, Bw};
}
dIR3 equilibrium_stellnaqs::del_contravariant(
    const IR3& position, double time) const {
  double r = position[IR3::u], theta = position[IR3::v];
  double phi = metric_->reduce_phi(position[IR3::w]);
  double coso = std::cos(theta), sino = std::sin(theta);
  double k = (*metric_->curvature())(phi);
  double k_prime = metric_->curvature()->derivative(phi);
  double d_u_Bw = 2.0*length_factor_*k*coso;
  double d_v_Bw = -2.0*length_factor_*k*r*sino;
  double d_w_Bw = 2.0*length_factor_*k_prime*r*coso;
  return {
      0.0, 0.0, 0.0, 
      iota_factor_*d_u_Bw, iota_factor_*d_v_Bw, iota_factor_*d_w_Bw,
      d_u_Bw, d_v_Bw, d_w_Bw};
}

}// end namespace gyronimo.
