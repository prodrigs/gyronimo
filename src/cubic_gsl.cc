// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @cubic_gsl.cc

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>

namespace gyronimo {

cubic_gsl::cubic_gsl(const dblock& x_range, const dblock& y_range)
    : spline_(nullptr), acc_(nullptr) {
  acc_ = gsl_interp_accel_alloc();
  spline_ = gsl_spline_alloc(gsl_interp_cspline, x_range.size());
  if (!spline_) error(
      __func__, __FILE__, __LINE__, "cannot allocate spline.", 1);
  gsl_spline_init(spline_, x_range.data(), y_range.data(), x_range.size());
}
cubic_gsl::~cubic_gsl() {
  if(spline_) gsl_spline_free(spline_);
  if(acc_) gsl_interp_accel_free(acc_);
}
double cubic_gsl::operator()(double x) const {
  return gsl_spline_eval(spline_, x, acc_);
}
double cubic_gsl::derivative(double x) const {
  return  gsl_spline_eval_deriv(spline_, x, acc_);
}
double cubic_gsl::derivative2(double x) const {
  return  gsl_spline_eval_deriv2(spline_, x, acc_);
}

} // end namespace gyronimo.
