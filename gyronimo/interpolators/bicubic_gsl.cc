// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @bicubic_gsl.cc, this file is part of ::gyronimo::

#include <algorithm>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/transpose.hh>
#include <gyronimo/interpolators/bicubic_gsl.hh>

namespace gyronimo {

bicubic_gsl::bicubic_gsl(
    const dblock& x_range, const dblock& y_range, const dblock& z_range,
    bool is_1st_faster, size_t periodic_size, size_t reflection_size)
    : spline_(nullptr), xacc_(nullptr), yacc_(nullptr) {
  if (periodic_size && reflection_size)  // if any, there can be only one!
      error(__func__, __FILE__, __LINE__,
          "incompatible boundary-condition request.", 1);
  if (periodic_size > y_range.size()/2 || reflection_size > y_range.size()/2)
      error(__func__, __FILE__, __LINE__,
          "boundary-condition extension is too large.", 1);
  xacc_ = gsl_interp_accel_alloc();
  yacc_ = gsl_interp_accel_alloc();
  if (!periodic_size && !reflection_size) {  // default, care only about order:
    spline_ = gsl_spline2d_alloc(
        gsl_interp2d_bicubic, x_range.size(), y_range.size());
    if (!spline_)
        error(__func__, __FILE__, __LINE__, "cannot allocate spline.", 1);
    gsl_spline2d_init(
        spline_, x_range.data(), y_range.data(),
        (is_1st_faster ?
            z_range.data() : transpose(z_range, y_range.size()).data()),
        x_range.size(), y_range.size());
  } else { // *Expensively* copies data into extended arrays:
    size_t nrows = periodic_size + reflection_size;  // ok, only one is not 0.
    const double* augmented_y =
        this->augment_y_range(y_range.data(), y_range.size(), nrows);
    const double* augmented_z = this->augment_z_range(
        (is_1st_faster ?
            z_range.data() : transpose(z_range, y_range.size()).data()),
        x_range.size(), y_range.size(),  // now, x is *always* faster.
        periodic_size, reflection_size);
    spline_ = gsl_spline2d_alloc(
        gsl_interp2d_bicubic, x_range.size(), y_range.size() + 2*nrows);
    if (!spline_) error(
        __func__, __FILE__, __LINE__, "cannot allocate spline.", 1);
    gsl_spline2d_init(spline_,
        x_range.data(), augmented_y, augmented_z,
            x_range.size(), y_range.size() + 2*nrows);
    delete [] augmented_z;
    delete [] augmented_y;
  }
}
bicubic_gsl::~bicubic_gsl() {
  if (spline_) gsl_spline2d_free(spline_);
  if (yacc_) gsl_interp_accel_free(yacc_);
  if (xacc_) gsl_interp_accel_free(xacc_);
}
double* bicubic_gsl::augment_z_range(
    const double* original, size_t nfast, size_t nslow,
    size_t periodic_size, size_t reflection_size) const {
  size_t nrows = periodic_size + reflection_size;  // ok, only one is not 0.
  size_t augmented_size = nfast*(nslow + 2*nrows);
  double* augmented = new double[augmented_size];
  if (!augmented) error(
      __func__, __FILE__, __LINE__, "cannot allocate augmented array.", 1);
  size_t offset = nrows*nfast;
  size_t original_size = nfast*nslow;
  std::copy(original, original + original_size, augmented + offset);
  if (periodic_size) {
    const double* lo = original + nfast;
    const double* hi = original + original_size - nfast;
    std::copy(hi - offset, hi, augmented);
    std::copy(lo, lo + offset, augmented + augmented_size - offset);
  } 
  if (reflection_size) {
    double* lo = augmented + offset;
    double* hi = augmented + augmented_size - offset - nfast;
    for (size_t row = 1; row <= reflection_size; row++) {
      std::copy(lo + row*nfast, lo + (row + 1)*nfast, lo - row*nfast);
      std::copy(hi - row*nfast, hi - (row - 1)*nfast, hi + row*nfast);
    }
  }
  return augmented;
}
double* bicubic_gsl::augment_y_range(
    const double* y_range, size_t ysize, size_t nrows) const {
  size_t augmented_size = ysize + 2*nrows;
  double* augmented = new double[augmented_size];
  if (!augmented) error(
      __func__, __FILE__, __LINE__, "cannot allocate augmented array.", 1);
  const double* bottom =  y_range + 1;
  const double* top = y_range + ysize - 1;
  std::copy(y_range, y_range + ysize, augmented + nrows);
  double period = y_range[ysize - 1] - y_range[0];
  std::transform(
      top - nrows, top, augmented, [period](double x){return x - period;});
  std::transform(
      bottom, bottom + nrows, augmented + augmented_size - nrows,
          [period](double x){return x + period;});
  return augmented;
}
double bicubic_gsl::operator()(double x, double y) const {
  return gsl_spline2d_eval(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial_u(double x, double y) const {
  return  gsl_spline2d_eval_deriv_x(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial_v(double x, double y) const {
  return gsl_spline2d_eval_deriv_y(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial2_uu(double x, double y) const {
  return  gsl_spline2d_eval_deriv_xx(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial2_uv(double x, double y) const {
  return  gsl_spline2d_eval_deriv_xy(spline_, x, y, xacc_, yacc_);
}
double bicubic_gsl::partial2_vv(double x, double y) const {
  return  gsl_spline2d_eval_deriv_yy(spline_, x, y, xacc_, yacc_);
}

} // end namespace gyronimo.
