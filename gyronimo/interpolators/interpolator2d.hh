// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @interpolator2d.hh

#ifndef GYRONIMO_INTERPOLATOR2D
#define GYRONIMO_INTERPOLATOR2D

#include <gyronimo/core/dblock.hh>

namespace gyronimo {

//! Access interface for 2d interpolators.
/*!
    Notice that this class only requires the **access** functionality to be
    implemented (evaluation, derivatives, second derivatives). Check the
    documentation of `interpolator1d` and `interpolator1d_factory` for details
    about the creation of specific interpolator objects by abstract code.
*/
class interpolator2d {
 public:
  virtual ~interpolator2d() {};
  virtual double operator()(double x, double y) const = 0;
  virtual double partial_u(double x, double y) const = 0;
  virtual double partial_v(double x, double y) const = 0;
  virtual double partial2_uu(double x, double y) const = 0;
  virtual double partial2_uv(double x, double y) const = 0;
  virtual double partial2_vv(double x, double y) const = 0;
};

//! Creation interface for 2d interpolators.
/*!
    Check the documentation of interpolator1d and interpolator1d_factory for
    details about the creation of specific interpolator objects by abstract
    code.
*/
class interpolator2d_factory {
 public:
  virtual interpolator2d* interpolate_data(
      const dblock& x_range,
      const dblock& y_range,
      const dblock& z_range) const = 0;
};

} // end namespace gyronimo.

#endif // GYRONIMO_INTERPOLATOR2D
