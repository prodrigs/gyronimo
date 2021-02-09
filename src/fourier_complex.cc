// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @fourier_complex.cc

#include <numeric>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/dblock.hh>
#include <gyronimo/interpolators/fourier_complex.hh>

namespace gyronimo {

fourier_complex::fourier_complex(
    const narray_t& u, const narray_t& dreal, const narray_t& dimag,
    const int mi, const int mf, const interpolator1d_factory* ifactory)
    : m_(mf - mi + 1),
      Areal_(mf - mi + 1, nullptr), Aimag_(mf - mi + 1, nullptr) {
  std::iota(std::begin(m_), std::end(m_), mi);
  this->build_interpolators(u, dreal, dimag, ifactory);
}
fourier_complex::fourier_complex(
    const narray_t& u, const narray_t& dreal, const narray_t& dimag,
    const narray_t& _m, const interpolator1d_factory* ifactory)
    : m_(std::begin(_m), std::end(_m)),
      Areal_(_m.size(), nullptr), Aimag_(_m.size(), nullptr) {
  this->build_interpolators(u, dreal, dimag, ifactory);
}
fourier_complex::~fourier_complex() {
  for (size_t p = 0; p < m_.size(); p++) {
    if (Areal_[p]) delete Areal_[p];
    if (Aimag_[p]) delete Aimag_[p];
  }
}
std::complex<double> fourier_complex::operator()(double u, double v) const {
  using namespace std::complex_literals;
  std::complex<double> sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++)
    sum += ((*Areal_[p])(u) + 1i*(*Aimag_[p])(u))*std::exp(1i*m_[p]*v);
  return sum;
}
std::complex<double> fourier_complex::partial_u(double u, double v) const {
  using namespace std::complex_literals;
  std::complex<double> sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++)
    sum += ((*Areal_[p]).derivative(u) +
        1i*(*Aimag_[p]).derivative(u))*std::exp(1i*m_[p]*v);
  return sum;
}
std::complex<double> fourier_complex::partial_v(double u, double v) const {
  using namespace std::complex_literals;
  std::complex<double> sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++)
    sum += 1i*m_[p]*((*Areal_[p])(u) + 1i*(*Aimag_[p])(u))*std::exp(1i*m_[p]*v);
  return sum;
}
std::complex<double> fourier_complex::partial2_uu(double u, double v) const {
  using namespace std::complex_literals;
  std::complex<double> sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++)
    sum += ((*Areal_[p]).derivative2(u) +
        1i*(*Aimag_[p]).derivative2(u))*std::exp(1i*m_[p]*v);
  return sum;
}
std::complex<double> fourier_complex::partial2_uv(double u, double v) const {
  using namespace std::complex_literals;
  std::complex<double> sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++)
    sum += 1i*m_[p]*((*Areal_[p]).derivative(u) +
        1i*(*Aimag_[p]).derivative(u))*std::exp(1i*m_[p]*v);
  return sum;
}
std::complex<double> fourier_complex::partial2_vv(double u, double v) const {
  using namespace std::complex_literals;
  std::complex<double> sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++)
    sum += -(m_[p]*m_[p])*(
      (*Areal_[p])(u) + 1i*(*Aimag_[p])(u))*std::exp(1i*m_[p]*v);
  return sum;
}
void fourier_complex::build_interpolators(
    const narray_t& u, const narray_t& dreal, const narray_t& dimag,
    const interpolator1d_factory* ifactory) {
  if (dreal.size() != dimag.size())
    error(__func__, __FILE__, __LINE__, "mismatched dreal and dimag.", 1);
  if (dreal.size() != u.size()*m_.size())
    error(__func__, __FILE__, __LINE__, "mismatched dreal, dimag, or u.", 1);
  for (size_t p = 0; p < m_.size(); p++) {
    std::slice index_range(p*u.size(), u.size(), 1);
    Areal_[p] = ifactory->interpolate_data(
        dblock_adapter(u), dblock_adapter<narray_t>(dreal[index_range]));
    Aimag_[p] = ifactory->interpolate_data(
        dblock_adapter(u), dblock_adapter<narray_t>(dimag[index_range]));
  }
}

} // end namespace gyronimo.
