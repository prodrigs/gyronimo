// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @equilibrium_circular.cc, this file is part of ::gyronimo::

#include <cmath>
#include <gyronimo/fields/equilibrium_circular.hh>

namespace gyronimo {

IR3 equilibrium_circular::contravariant(
    const IR3& position, double time) const {
  double m = this->magnitude(position, time);
  IR3 A = this->contravariant_versor(position, time);
  return {m * A[IR3::u], m * A[IR3::v], m * A[IR3::w]};
}
dIR3 equilibrium_circular::del_contravariant(
    const IR3& position, double time) const {
  double R0 = metric_->major_radius();
  double eps = metric_->iaspect_ratio();
  double r = position[IR3::u], theta = position[IR3::v];
  double q = q_(r), qprime = qprime_(r);
  double dRdr = eps*std::cos(theta);
  double dRdtheta = -r*eps*std::sin(theta);
  double R = 1.0 + r*dRdr;
  double dB_vu = -(qprime*R + q*dRdr)/(R0*q*q*R*R);
  double dB_vv = -q*dRdtheta/(R0*q*q*R*R);
  double dB_wu = -2.0*dRdr/(R0*R*R*R);
  double dB_wv = -2.0*dRdtheta/(R0*R*R*R);
  return {0.0, 0.0, 0.0, dB_vu, dB_vv, 0.0, dB_wu, dB_wv, 0.0};
}
IR3 equilibrium_circular::covariant(const IR3& position, double time) const {
  double m = this->magnitude(position, time);
  IR3 A = this->covariant_versor(position, time);
  return {m*A[IR3::u], m*A[IR3::v], m*A[IR3::w]};
}
double equilibrium_circular::magnitude(const IR3& position, double time) const {
  double eps_r = metric_->iaspect_ratio() * position[IR3::u];
  double q = q_(position[IR3::u]);
  double l = std::sqrt(q*q + eps_r*eps_r);
  double R = 1.0 + eps_r*std::cos(position[IR3::v]);
  return l/(q*R);
}
IR3 equilibrium_circular::contravariant_versor(
    const IR3& position, double time) const {
  double R0 = metric_->major_radius();
  double r = position[IR3::u], theta = position[IR3::v];
  double eps_r = metric_->iaspect_ratio()*r;
  double q = q_(r);
  double R = 1.0 + eps_r*std::cos(theta);
  double aux = 1.0/(R0*R*std::sqrt(q*q + eps_r*eps_r));
  return {0.0, R*aux, q*aux};
}
IR3 equilibrium_circular::covariant_versor(
    const IR3& position, double time) const {
  IR3 b = this->contravariant_versor(position, time);
  return this->metric()->to_covariant(b, position);
}
IR3 equilibrium_circular::curl(const IR3& position, double time) const {
  double J = this->metric()->jacobian(position);
  dIR3 dB = this->del_covariant(position, time);
  return {0.0, 0.0, dB[dIR3::vu]/J};
}
IR3 equilibrium_circular::del_magnitude(
    const IR3& position, double time) const {
  double r = position[IR3::u], theta = position[IR3::v];
  double eps = metric_->iaspect_ratio();
  double eps_r = eps*r;
  double q = q_(r);
  double qprime = qprime_(r);
  double l = std::sqrt(q*q + eps_r*eps_r);
  double lprime = (q*qprime + eps*eps_r)/l;
  double duR = eps*std::cos(theta);
  double R = 1.0 + r*duR;
  double aux = std::pow(q*R, -2);
  return {
      (q*R*lprime - l*(R*qprime + q*duR))*aux,
       q*l*eps_r*std::sin(theta)*aux, 0.0};
}

} // end namespace gyronimo.
