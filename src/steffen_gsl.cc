// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @steffen_gsl.cc

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/steffen_gsl.hh>

namespace gyronimo {

steffen_gsl::steffen_gsl(const dblock& x_range, const dblock& y_range)
    : spline_(nullptr), acc_(nullptr) {
  acc_ = gsl_interp_accel_alloc();
  spline_ = gsl_spline_alloc(gsl_interp_steffen, x_range.size());
  if (!spline_) error(
      __func__, __FILE__, __LINE__, "cannot allocate spline.", 1);
  gsl_spline_init(spline_, x_range.data(), y_range.data(), x_range.size());
}
steffen_gsl::~steffen_gsl() {
  if(spline_) gsl_spline_free(spline_);
  if(acc_) gsl_interp_accel_free(acc_);
}
double steffen_gsl::operator()(double x) const {
  return gsl_spline_eval(spline_, x, acc_);
}
double steffen_gsl::derivative(double x) const {
  return  gsl_spline_eval_deriv(spline_, x, acc_);
}
double steffen_gsl::derivative2(double x) const {
  return  gsl_spline_eval_deriv2(spline_, x, acc_);
}

} // end namespace gyronimo.
