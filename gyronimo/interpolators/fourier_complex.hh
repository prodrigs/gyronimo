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

// @fourier_complex.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_FOURIER_COMPLEX
#define GYRONIMO_FOURIER_COMPLEX

#include <vector>
#include <complex>
#include <valarray>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo {

//! Fourier representation of a complex field over the plane.
/*!
    Implements the formula @f$ f(u,v) = \sum_{m=m_i}^{m_f} f_m(u) e^{i m v} @f$,
    where @f$ f_m(u) @f$ are complex-valued functions of the real argument u.
    These functions are to be built by interpolating data supplied to the
    constructor. In the constructors, `u` is an array of grid locations, while
    `dreal` and `dimag` are the corresponding real and imaginary samples of the
    field for each harmonic. If `N=n.size()` and `M=(mf-mi+1)=_m.size()`, then
    one must have `dreal.size()=dimag.size()=N*M`. The specific type of
    interpolator to use is controlled by the object `ifactory` provided to the
    constructors. Although sharing a rather similar interface, `fourier_complex`
    is **not** derived from `interpolator2d` because its members return complex
    numbers instead of real ones.
    
    @todo Change constructor input from narray_t to general containers.
*/
class fourier_complex {
 public:
  typedef std::valarray<double> narray_t;
  fourier_complex(
      const narray_t& u, const narray_t& dreal, const narray_t& dimag,
      const int mi, const int mf, const interpolator1d_factory* ifactory);
  fourier_complex(
      const narray_t& u, const narray_t& dreal, const narray_t& dimag,
      const narray_t& _m, const interpolator1d_factory* ifactory);
  ~fourier_complex();
  std::complex<double> operator()(double u, double v) const;
  std::complex<double> partial_u(double u, double v) const;
  std::complex<double> partial_v(double u, double v) const;
  std::complex<double> partial2_uu(double u, double v) const;
  std::complex<double> partial2_uv(double u, double v) const;
  std::complex<double> partial2_vv(double u, double v) const;

 private:
  std::vector<double> m_;
  std::vector<interpolator1d*> Areal_;
  std::vector<interpolator1d*> Aimag_;
  void build_interpolators(
      const narray_t& u, const narray_t& dreal, const narray_t& dimag,
      const interpolator1d_factory* ifactory);
};

} // end namespace gyronimo.

#endif // GYRONIMO_FOURIER_COMPLEX
