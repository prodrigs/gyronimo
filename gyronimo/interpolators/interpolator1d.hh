// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @interpolator1d.hh

#ifndef GYRONIMO_INTERPOLATOR1D
#define GYRONIMO_INTERPOLATOR1D

#include <gyronimo/core/dblock.hh>

namespace gyronimo {

//! Access interface for 1d interpolators.
/*!
    Notice that this class only requires the **access** functionality to be
    implemented (evaluation, derivative, second derivative). The **creation** of
    specific interpolator objects (corresponding to classes derived from
    interpolator1d) by abstract code is to be handled by classes derived from
    interpolator1d_factory, whose documentation should be checked for more
    details.
*/
class interpolator1d {
 public:
  virtual ~interpolator1d() {};
  virtual double operator()(double x) const = 0;
  virtual double derivative(double x) const = 0;
  virtual double derivative2(double x) const = 0;
};

//! Creation interface for 1d interpolators.
/*!
    Sometimes, a piece of general and abstract code that only knows about the
    interface of `interpolator1d` needs to **create** a specific interpolator
    object of a specific kind (i.e., a class derived from `interpolator1d`)
    requested by the user. This can be achieved using a pointer to
    `interpolator1d_factory` pointing to an object that overrides the member
    interpolate_data with the actual creation of the intended interpolator. To
    benefit from this level of abstraction, every interpolator class derived
    from interpolator1d **must** have a corresponding class derived from
    `interpolator1d_factory`.
*/
class interpolator1d_factory {
 public:
  virtual interpolator1d* interpolate_data(
      const dblock& x_range, const dblock& y_range) const = 0;
};

} // end namespace gyronimo.

#endif // GYRONIMO_INTERPOLATOR1D
