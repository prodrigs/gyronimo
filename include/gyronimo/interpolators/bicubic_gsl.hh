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
/*!
    The grid samples in `x_range` and `y_range` may be non-uniform (there are no
    performance gains by letting them be uniform), but `x_range.size() *
    y_range.size() = z_range.size()`. Which variable changes faster in `z_range`
    is described by the boolean argument `is_1st_faster`. By default, natural
    boundary conditions are assumed for both dimensions. However, support is
    provided for periodic and reflection boundary conditions for the **second
    variable only** by extending the domain by `periodic_size` or
    `reflection_size` samples on each side.
 */
class bicubic_gsl : public interpolator2d {
 public:
  bicubic_gsl(
      const dblock& x_range, const dblock& y_range, const dblock& z_range,
      bool is_1st_faster, size_t periodic_size = 0, size_t reflection_size = 0);
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
  double* augment_y_range(
      const double* y_range, size_t ysize, size_t nrows) const;
  double* augment_z_range(
      const double* original, size_t nfast, size_t nslow,
      size_t periodic_size, size_t reflection_size) const;
};

class bicubic_gsl_factory : public interpolator2d_factory {
 public:
  bicubic_gsl_factory(
      bool is_1st_faster, size_t periodic_size = 0, size_t reflection_size = 0)
      : is_1st_faster_(is_1st_faster),
        periodic_size_(periodic_size), reflection_size_(reflection_size) {};
  virtual interpolator2d* interpolate_data(
      const dblock& x_range,
      const dblock& y_range,
      const dblock& z_range) const override {
    return new bicubic_gsl(x_range, y_range, z_range,
        is_1st_faster_, periodic_size_, reflection_size_);
  };
 private:
  bool is_1st_faster_;
  size_t periodic_size_, reflection_size_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_BICUBIC_GSL
