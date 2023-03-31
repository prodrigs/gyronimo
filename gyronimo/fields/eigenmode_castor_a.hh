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

// @eigenmode_castor_a.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_EIGENMODE_CASTOR_A
#define GYRONIMO_EIGENMODE_CASTOR_A

#include <gyronimo/fields/eigenmode_castor_b.hh>

namespace gyronimo {

//! Vector-potential eigenvector from a `CASTOR` output file.
/*!
    The vector-potential, divided by `m_factor` (SI), is the adimensional
    perturbation that `CASTOR` relates with the magnetic-field via
    @f$\mathbf{B}_\mathrm{cas} = \tilde{\nabla} \times
    \mathbf{A}_\mathrm{cas}@f$, i.e., the adimensional equation 2.18 in [W.
    Kerner *et al*., J. Comput. Phys. **142**, 271 (1998)], where @f$
    \tilde{\nabla} = R_0 \nabla @f$, and @f$R_0@f$ is the magnetic axis radius.
    Its covariant components are represented as:
    @f[
      A_k^\mathrm{cas}(s, \chi, \phi, t) = e^{\lambda t + i n \phi}
          \sum_m \hat{A}_{k,m}^\mathrm{cas}(s) e^{i m \chi}.
    @f]
    Each @f$ \hat{A}^\mathrm{cas}_{k,m}(s) @f$ is a `fourier_complex` object
    with underlying `interpolator1d` type set by `ifactory`. The time
    normalisation is set to the Alfven time @f$R_0/v_A^0@f$, with @f$v_A^0@f$
    the Alfven velocity on axis that is supplied to the constructor via
    `v_alfven` (SI). At construction, the field @f$\mathbf{A}_\mathrm{cas}@f$ is
    multiplied by a given constant (returned by `native_factor()`) devised in
    order to set its maximum magnitude over the poloidal cross section (i.e.,
    @f$\phi = 0, t = 0@f$) to unity. An alternative constructor allows the
    normalisation factors to be automatically deduced in order to enforce
    compatibility with a parent `eigenmode_castor_b` object.

    @todo higher-order derivatives from the interpolators seem to produce noisy
    magnetic fields; implement derivatives produced by castor itself (from the
    internal hermite elements, not yet available in parser_castor).

    @todo move the normalisation done in the constructor into a code block
    common with all elements of the `eigenmode_castor_x` family. Extend
    automatic normalisation construction to other eigenmode types.
*/
class eigenmode_castor_a : public IR3field {
 public:
  eigenmode_castor_a(
      double m_factor, double v_alfven,
      const parser_castor *p, const metric_helena *g,
      const interpolator1d_factory* ifactory);
  eigenmode_castor_a(
      const eigenmode_castor_b* parent, const interpolator1d_factory* ifactory);
  virtual ~eigenmode_castor_a() override {};

  virtual IR3 covariant(const IR3& position, double time) const override;
  virtual IR3 contravariant(const IR3& position, double time) const override;

  const parser_castor* parser() const {return parser_;};
  double native_factor() const {return native_factor_;};
  double v_alfven() const {return metric_->parser()->rmag()/this->t_factor();};

 private:
  double native_factor_;
  const parser_castor *parser_;
  const metric_helena *metric_;
  const std::complex<double> eigenvalue_, i_n_tor_;
  fourier_complex tildeA1_, tildeA2_, tildeA3_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_EIGENMODE_CASTOR_A
