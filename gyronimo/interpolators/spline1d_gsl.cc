// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @spline1d_gsl.cc

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/spline1d_gsl.hh>

namespace gyronimo {

spline1d_gsl::spline1d_gsl() : spline_(nullptr), acc_(nullptr) {
  acc_ = gsl_interp_accel_alloc();
}
spline1d_gsl::~spline1d_gsl() {
  if (acc_) delete(acc_);
}
double spline1d_gsl::operator()(double x) const {
  return gsl_spline_eval(spline_, x, acc_);
}
double spline1d_gsl::derivative(double x) const {
  return  gsl_spline_eval_deriv(spline_, x, acc_);
}
double spline1d_gsl::derivative2(double x) const {
  return  gsl_spline_eval_deriv2(spline_, x, acc_);
}

} // end namespace gyronimo.
