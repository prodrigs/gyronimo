// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues.

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

// @eigenmode_castor_b.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_EIGENMODE_CASTOR_B
#define GYRONIMO_EIGENMODE_CASTOR_B

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/parsers/parser_castor.hh>
#include <gyronimo/metrics/metric_helena.hh>
#include <gyronimo/interpolators/interpolator1d.hh>
#include <gyronimo/interpolators/fourier_complex.hh>

namespace gyronimo {

//! Magnetic-field eigenvector from a `CASTOR` output file.
/*!
    The magnetic-field, divided by `m_factor` (SI), is the adimensional magnetic
    perturbation as defined by `CASTOR` from the vector-potential via
    @f$\mathbf{B}_\mathrm{cas} = \tilde{\nabla} \times
    \mathbf{A}_\mathrm{cas}@f$, i.e., the adimensional equation 2.18 in [W.
    Kerner *et al*., J. Comput. Phys. **142**, 271 (1998)], where @f$
    \tilde{\nabla} = R_0 \nabla @f$, and @f$R_0@f$ is the magnetic axis radius.
    On the other hand, covariant components of the adimensional vector-potential
    as defined by `CASTOR` are represented as:
    @f[
      A_k^\mathrm{cas}(s, \chi, \phi, t) = e^{\lambda t + i n \phi}
          \sum_m \hat{A}_{k,m}^\mathrm{cas}(s) e^{i m \chi}.
    @f]
    Each @f$ \hat{A}^\mathrm{cas}_{k,m}(s) @f$ is a `fourier_complex` object
    with underlying `interpolator1d` type set by `ifactory`. The time
    normalisation is set to the Alfven time @f$R_0/v_A^0@f$, with @f$v_A^0@f$
    the Alfven velocity on axis that is supplied to the constructor via
    `v_alfven` (SI). At construction, the field @f$\mathbf{B}_\mathrm{cas}@f$ is
    multiplied by a given constant (returned by `native_factor()`) devised in
    order to set its maximum magnitude over the poloidal cross section (i.e.,
    @f$\phi = 0, t = 0@f$) to unity.

    @todo higher-order derivatives from the interpolators seem to produce noisy
    magnetic fields; implement derivatives produced by castor itself (from the
    internal hermite elements, not yet available in parser_castor).

    @todo move the normalisation done in the constructor into a code block
    common with all elements of the `eigenmode_castor_x` family, along with the
    common interface.
*/
class eigenmode_castor_b : public IR3field_c1 {
 public:
  eigenmode_castor_b(
      double m_factor, double v_alfven,
      const parser_castor *p, const metric_helena *g,
      const interpolator1d_factory* ifactory);
  virtual ~eigenmode_castor_b() override {};

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override;

  const parser_castor* parser() const {return parser_;};
  double native_factor() const {return native_factor_;};
  double v_alfven() const {return metric_->parser()->rmag()/this->t_factor();};

 private:
  double native_factor_;
  const parser_castor *parser_;
  const metric_helena *metric_;
  double n_tor_squared_;
  std::complex<double> eigenvalue_, i_n_tor_;
  fourier_complex tildeA1_, tildeA2_, tildeA3_;

  inline std::complex<double> exp_wt_nphi(double time, double phi) const;
  inline std::complex<double> d1A2(double s, double chi) const;
  inline std::complex<double> d1A3(double s, double chi) const;
  inline std::complex<double> d2A1(double s, double chi) const;
  inline std::complex<double> d2A3(double s, double chi) const;
  inline std::complex<double> d3A1(double s, double chi) const;
  inline std::complex<double> d3A2(double s, double chi) const;
  inline std::complex<double> d11A2(double s, double chi) const;
  inline std::complex<double> d11A3(double s, double chi) const;
  inline std::complex<double> d21A1(double s, double chi) const;
  inline std::complex<double> d21A3(double s, double chi) const;
  inline std::complex<double> d31A1(double s, double chi) const;
  inline std::complex<double> d31A2(double s, double chi) const;
  inline std::complex<double> d12A2(double s, double chi) const;
  inline std::complex<double> d12A3(double s, double chi) const;
  inline std::complex<double> d22A1(double s, double chi) const;
  inline std::complex<double> d22A3(double s, double chi) const;
  inline std::complex<double> d32A1(double s, double chi) const;
  inline std::complex<double> d32A2(double s, double chi) const;
  inline std::complex<double> d13A2(double s, double chi) const;
  inline std::complex<double> d13A3(double s, double chi) const;
  inline std::complex<double> d23A1(double s, double chi) const;
  inline std::complex<double> d23A3(double s, double chi) const;
  inline std::complex<double> d33A1(double s, double chi) const;
  inline std::complex<double> d33A2(double s, double chi) const;
};

} // end namespace gyronimo.

#endif // GYRONIMO_EIGENMODE_CASTOR_B
