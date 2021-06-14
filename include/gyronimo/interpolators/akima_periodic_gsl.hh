// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @akima_periodic_gsl.hh

#ifndef GYRONIMO_AKIMA_PERIODIC_GSL
#define GYRONIMO_AKIMA_PERIODIC_GSL

#include <gsl/gsl_spline.h>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo {

//! Periodic Akima spline by [GSL](https://www.gnu.org/software/gsl).
class akima_periodic_gsl : public interpolator1d {
 public:
  akima_periodic_gsl(const dblock& x_range, const dblock& y_range);
  virtual ~akima_periodic_gsl() override;

  double operator()(double x) const override;
  double derivative(double x) const override;
  double derivative2(double x) const override;

 private:
  gsl_spline *spline_;
  gsl_interp_accel *acc_;
};

class akima_periodic_gsl_factory : public interpolator1d_factory {
 public:
  virtual interpolator1d* interpolate_data(
      const dblock& x_range, const dblock& y_range) const override {
    return new akima_periodic_gsl(x_range, y_range);
  };
};

} // end namespace gyronimo.

#endif // GYRONIMO_AKIMA_PERIODIC_GSL
