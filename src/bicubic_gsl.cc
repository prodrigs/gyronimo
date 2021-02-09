// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @bicubic_gsl.cc

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/bicubic_gsl.hh>

namespace gyronimo {

bicubic_gsl::bicubic_gsl(
    const dblock& x_range, const dblock& y_range, const dblock& z_range)
    : spline_(nullptr), xacc_(nullptr), yacc_(nullptr) {
  xacc_ = gsl_interp_accel_alloc();
  yacc_ = gsl_interp_accel_alloc();
  spline_ = gsl_spline2d_alloc(
      gsl_interp2d_bicubic, x_range.size(), y_range.size());
  if (!spline_) error(
      __func__, __FILE__, __LINE__, "cannot allocate spline.", 1);
  gsl_spline2d_init(
      spline_,
      x_range.data(), y_range.data(), z_range.data(),
      x_range.size(), y_range.size());
}
bicubic_gsl::~bicubic_gsl() {
  if (spline_) gsl_spline2d_free(spline_);
  if (yacc_) gsl_interp_accel_free(yacc_);
  if (xacc_) gsl_interp_accel_free(xacc_);
}
double bicubic_gsl::operator()(double x, double y) const {
  return gsl_spline2d_eval(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial_u(double x, double y) const {
  return  gsl_spline2d_eval_deriv_x(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial_v(double x, double y) const {
  return gsl_spline2d_eval_deriv_y(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial2_uu(double x, double y) const {
  return  gsl_spline2d_eval_deriv_xx(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial2_uv(double x, double y) const {
  return  gsl_spline2d_eval_deriv_xy(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial2_vv(double x, double y) const {
  return  gsl_spline2d_eval_deriv_yy(spline_, x, y, xacc_, yacc_);
}

} // end namespace gyronimo.
