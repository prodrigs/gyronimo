// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @spline1d_gsl.hh

#ifndef GYRONIMO_SPLINE1D_GSL
#define GYRONIMO_SPLINE1D_GSL

#include <gsl/gsl_spline.h>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo {

//! Base class for 1d splines by [GSL](https://www.gnu.org/software/gsl).
class spline1d_gsl : public interpolator1d {
 public:
  spline1d_gsl();
  virtual ~spline1d_gsl() override;

  double operator()(double x) const final;
  double derivative(double x) const final;
  double derivative2(double x) const final;

 protected:
  gsl_spline *spline_;
  gsl_interp_accel *acc_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_SPLINE1D_GSL
