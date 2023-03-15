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

// @eigenmode_castor_e.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_EIGENMODE_CASTOR_E
#define GYRONIMO_EIGENMODE_CASTOR_E

#include <gyronimo/fields/eigenmode_castor_a.hh>

namespace gyronimo {

//! Electric-field eigenvector from `CASTOR` ceig output file.
/*!
    The electric-field **covariant** components are extracted as @f$\mathbf{E} =
    - \partial_t \mathbf{A}@f$ from the vector potential stored in a
    `eigenmode_castor_a` object. The time normalisation `t_factor` **must** be
    set to the ratio @f$R_0/v_A(0)@f$ of the axis radius to the on-axis Alfven
    velocity. The field is normalised internally by setting its maximum
    magnitude over the poloidal cross section (@f$\phi = 0@f$) to the value
    `m_factor` (SI).
*/

class eigenmode_castor_e : public IR3field {
 public:
  eigenmode_castor_e(
      double m_factor, double t_factor,
      const parser_castor *p, const metric_helena *g,
      const interpolator1d_factory* ifactory);
  eigenmode_castor_e(
      double m_factor, double t_factor,
      const parser_castor *p, const metric_helena *g,
      const interpolator1d_factory* ifactory, 
      double norm_factor);
  virtual ~eigenmode_castor_e() override {};

  virtual IR3 covariant(const IR3& position, double time) const override;
  virtual IR3 contravariant(const IR3& position, double time) const override;

  const parser_castor* parser() const {return parser_;};
  double get_norm_factor(){return norm_factor_;};

 private:
  double norm_factor_;
  const parser_castor *parser_;
  const metric_helena *metric_;
  std::complex<double> eigenvalue_, i_n_tor_;
  fourier_complex tildeA1_, tildeA2_, tildeA3_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_EIGENMODE_CASTOR_E
