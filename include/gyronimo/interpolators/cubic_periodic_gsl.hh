// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @cubic_periodic_gsl.hh

#ifndef GYRONIMO_CUBIC_PERIODIC_GSL
#define GYRONIMO_CUBIC_PERIODIC_GSL

#include <gyronimo/interpolators/spline1d_gsl.hh>

namespace gyronimo {

//! Periodic cubic spline by [GSL](https://www.gnu.org/software/gsl).
class cubic_periodic_gsl : public spline1d_gsl {
 public:
  cubic_periodic_gsl(const dblock& x_range, const dblock& y_range);
  virtual ~cubic_periodic_gsl() final;
};

class cubic_periodic_gsl_factory : public interpolator1d_factory {
 public:
  virtual interpolator1d* interpolate_data(
      const dblock& x_range, const dblock& y_range) const final {
    return new cubic_periodic_gsl(x_range, y_range);
  };
};

} // end namespace gyronimo.

#endif // GYRONIMO_CUBIC_PERIODIC_GSL
