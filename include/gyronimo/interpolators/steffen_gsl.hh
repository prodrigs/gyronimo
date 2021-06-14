// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @steffen_gsl.hh

#ifndef GYRONIMO_STEFFEN_GSL
#define GYRONIMO_STEFFEN_GSL

#include <gyronimo/interpolators/spline1d_gsl.hh>

namespace gyronimo {

//! Natural Steffen spline by [GSL](https://www.gnu.org/software/gsl).
class steffen_gsl : public spline1d_gsl {
 public:
  steffen_gsl(const dblock& x_range, const dblock& y_range);
  virtual ~steffen_gsl() final;
};

class steffen_gsl_factory : public interpolator1d_factory {
 public:
  virtual interpolator1d* interpolate_data(
      const dblock& x_range, const dblock& y_range) const final {
    return new steffen_gsl(x_range, y_range);
  };
};

} // end namespace gyronimo.

#endif // GYRONIMO_STEFFEN_GSL
