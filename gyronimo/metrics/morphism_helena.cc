// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues.

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

#include <gyronimo/metrics/morphism_helena.hh>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gyronimo/core/multiroot.hh>

namespace gyronimo{

morphism_helena::morphism_helena(
    const parser_helena *p, const interpolator2d_factory *ifactory)
    : parser_(p), R_(nullptr), z_(nullptr), 
	R0_(p->rmag()), squaredR0_(p->rmag()*p->rmag()) {
  double Rgeo = p->rgeo();
  double a = p->eps()*Rgeo;
  dblock_adapter s_range(p->s()), chi_range(p->chi());
  R_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(
          parser_helena::narray_type(a*p->x() + Rgeo)));
  z_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(
          parser_helena::narray_type(a*p->y() + 0.0)));

  guu_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g11()));
  guv_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g12()));
  gvv_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g22()));
  gww_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g33()));
}
morphism_helena::~morphism_helena() {
  if(R_) delete R_;
  if(z_) delete z_;
  if(guu_) delete guu_;
  if(guv_) delete guv_;
  if(gvv_) delete gvv_;
  if(gww_) delete gww_;
}
IR3 morphism_helena::operator()(const IR3& q) const {
  double s = q[IR3::u], chi = reduce_2pi(q[IR3::v]), phi = q[IR3::w];
  double R = (*R_)(s, chi);
  return {R*std::cos(phi), -R*std::sin(phi), (*z_)(s, chi)};
}
// int morphism_helena::root_f(const gsl_vector *q, void *target, gsl_vector *f) {
// //   double u = reduce(gsl_vector_get(q, 0), 1);
//   reduce2(q->data[0], q->data[1]);
//   double u = q->data[0];
//   double v = reduce(gsl_vector_get(q, 1), 2*std::numbers::pi);
//   double target_R = ((struct root_target*) target)->R;
//   double target_z = ((struct root_target*) target)->z;
//   interpolator2d* interpolator_R = ((struct root_target*) target)->i_R;
//   interpolator2d* interpolator_z = ((struct root_target*) target)->i_z;
//   gsl_vector_set(f, 0, target_R - (*interpolator_R)(u, v));
//   gsl_vector_set(f, 1, target_z - (*interpolator_z)(u, v));
//   return GSL_SUCCESS;
// }
IR3 morphism_helena::inverse(const IR3& X) const {
  typedef std::array<double, 2> IR2;
  double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
  double R = std::sqrt(x*x + y*y);
  std::function<IR2(const IR2&)> zero_function =
      [&](const IR2& args) {
        auto [s, chi] = reflection_past_axis(args[0], args[1]);
        return IR2({(*R_)(s, chi) - R, (*z_)(s, chi) - z});
      };
  IR2 guess = {0.5, std::atan2(z, R - parser_->rmag())};
  IR2 roots = multiroot(1.0e-15, 100)(zero_function, guess);
  auto [s, chi] = reflection_past_axis(roots[0], roots[1]);
  return {s, chi, std::atan2(-y, x)};
}
dIR3 morphism_helena::del(const IR3& q) const {
  double s = q[IR3::u], chi = reduce_2pi(q[IR3::v]), phi = q[IR3::w];
  double R = (*R_)(s, chi),
      Ru = (*R_).partial_u(s, chi), Rv = (*R_).partial_v(s, chi);
  double cos = std::cos(phi), sin = std::sin(phi);
  return {Ru*cos, Rv*cos, -R*sin, -Ru*sin, -Rv*sin, -R*cos,
      (*z_).partial_u(s, chi), (*z_).partial_v(s, chi), 0.0};
}
double morphism_helena::jacobian(const IR3 &q) const {
	double s = q[IR3::u], chi = reduce_2pi(q[IR3::v]);
	double guu = squaredR0_ * (*guu_)(s, chi);
	double guv = squaredR0_ * (*guv_)(s, chi);
	double gvv = squaredR0_ * (*gvv_)(s, chi);
	double gww = squaredR0_ * (*gww_)(s, chi);
	return std::sqrt((guu*gvv - guv*guv)*gww);
}
dSM3 morphism_helena::g_del(const IR3& q) const {
	double s = q[IR3::u], chi = reduce_2pi(q[IR3::v]);
	return {
		squaredR0_*(*guu_).partial_u(s, chi),
		squaredR0_*(*guu_).partial_v(s, chi), 0.0, // d_i g_uu
		squaredR0_*(*guv_).partial_u(s, chi),
		squaredR0_*(*guv_).partial_v(s, chi), 0.0, // d_i g_uv
		0.0, 0.0, 0.0, // d_i g_uw
		squaredR0_*(*gvv_).partial_u(s, chi),
		squaredR0_*(*gvv_).partial_v(s, chi), 0.0, //d_i g_vv
		0.0, 0.0, 0.0, // d_i g_vw
		squaredR0_*(*gww_).partial_u(s, chi),
		squaredR0_*(*gww_).partial_v(s, chi), 0.0}; // d_i g_ww
}

// //! Reduces an arbitrary angle chi to the interval [0:pi].
// double morphism_helena::reduce_chi(double chi) const {
//   double rchi = reduce(chi, 2*std::numbers::pi);
//   if(parser_->is_symmetric() && rchi > std::numbers::pi)
//       rchi = 2*std::numbers::pi - rchi;
//   return rchi;
// }

std::tuple<double, double> morphism_helena::reflection_past_axis(
    double s, double chi) {
  if(s < 0)
    return {-s, reduce_2pi(chi + std::numbers::pi)};
  else
    return {s, chi};
}
double morphism_helena::reduce_2pi(double x) {
  double l = 2*std::numbers::pi;
  return (x -= l*std::floor(x/l));
}

} // end namespace gyronimo
