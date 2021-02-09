// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @bicubic_gsl.hh

#ifndef GYRONIMO_BICUBIC_GSL
#define GYRONIMO_BICUBIC_GSL

#include <gsl/gsl_spline2d.h>
#include <gyronimo/interpolators/interpolator2d.hh>

namespace gyronimo {

//! Bicubic spline using the [GSL library](https://www.gnu.org/software/gsl).
class bicubic_gsl : public interpolator2d {
 public:
  bicubic_gsl(
      const dblock& x_range, const dblock& y_range, const dblock& z_range);
  virtual ~bicubic_gsl();

  virtual double operator()(double x, double y) const override;
  virtual double partial_u(double x, double y) const override;
  virtual double partial_v(double x, double y) const override;
  virtual double partial2_uu(double x, double y) const override;
  virtual double partial2_uv(double x, double y) const override;
  virtual double partial2_vv(double x, double y) const override;

 private:
  gsl_spline2d *spline_;
  gsl_interp_accel *xacc_;
  gsl_interp_accel *yacc_;
};

class bicubic_gsl_factory : public interpolator2d_factory {
 public:
  virtual interpolator2d* interpolate_data(
      const dblock& x_range,
      const dblock& y_range,
      const dblock& z_range) const override {
    return new bicubic_gsl(x_range, y_range, z_range);
  };
};

} // end namespace gyronimo.

#endif // GYRONIMO_BICUBIC_GSL
