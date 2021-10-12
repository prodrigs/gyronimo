// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @akima_gsl.hh

#ifndef GYRONIMO_AKIMA_GSL
#define GYRONIMO_AKIMA_GSL

#include <gyronimo/interpolators/spline1d_gsl.hh>

namespace gyronimo {

//! Natural Akima spline by [GSL](https://www.gnu.org/software/gsl).
class akima_gsl : public spline1d_gsl {
 public:
  akima_gsl(const dblock& x_range, const dblock& y_range);
  virtual ~akima_gsl() final;
};

class akima_gsl_factory : public interpolator1d_factory {
 public:
  virtual interpolator1d* interpolate_data(
      const dblock& x_range, const dblock& y_range) const final {
    return new akima_gsl(x_range, y_range);
  };
};

} // end namespace gyronimo.

#endif // GYRONIMO_AKIMA_GSL
