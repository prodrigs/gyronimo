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

// @eigenmode_castor_a.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_EIGENMODE_CASTOR_A
#define GYRONIMO_EIGENMODE_CASTOR_A

#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/parsers/parser_castor.hh>
#include <gyronimo/metrics/metric_helena.hh>
#include <gyronimo/interpolators/interpolator1d.hh>
#include <gyronimo/interpolators/fourier_complex.hh>

namespace gyronimo {

//! Vector-potential eigenvector from `CASTOR` ceig output file.
/*!
    The **covariant** components of the potential vector are represented as:
    @f[
      A_k(s, \chi, \phi, t) =
          e^{i \omega t + i n \phi} \sum_m \tilde{A}_k^m(s) e^{i m \chi}.
    @f]
    Each @f$ \tilde{A}_k^m(s) @f$ is a `fourier_complex` object with underlying
    `interpolator1d` type set by `ifactory`. The time normalisation `t_factor`
    is the ratio @f$v_A(0)/R_0@f$ of the on-axis Alfven velocity to the axis
    radius. The field is normalised to its maximum magnitude over the low-field
    side (@f$\chi = 0@f$), its actual amplitude being set by `m_factor` (SI).
*/
class eigenmode_castor_a : public IR3field {
 public:
  eigenmode_castor_a(
      double m_factor, double t_factor,
      const parser_castor *p, const metric_helena *g,
      const interpolator1d_factory* ifactory);
  virtual ~eigenmode_castor_a() override {};

  virtual IR3 covariant(const IR3& position, double time) const override;
  virtual IR3 contravariant(const IR3& position, double time) const override;

  const parser_castor* parser() const {return parser_;};

 private:
  double norm_factor_;
  const parser_castor *parser_;
  const metric_helena *metric_;
  std::complex<double> w_, i_n_tor_;
  fourier_complex tildeA1_, tildeA2_, tildeA3_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_EIGENMODE_CASTOR_A
