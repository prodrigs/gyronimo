// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022, 2023 Jorge Ferreira and Paulo Rodrigues.

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

// @metric_vmec.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_VMEC
#define GYRONIMO_METRIC_VMEC

#include <memory>
#include <complex>
#include <gyronimo/interpolators/interpolator1d.hh>
#include <gyronimo/metrics/metric_covariant.hh>
#include <gyronimo/parsers/parser_vmec.hh>

namespace gyronimo {

//! Covariant metric in `VMEC` curvilinear coordinates.
/*!
    Builds the metric from the information provided by a `parser_vmec` object.
    The right-handed coordinates are the poloidal flux normalised to its
    boundary value (`u`, or `VMEC` @f$\chi@f$), the toroidal angle (`v`, or
    `VMEC` @f$\zeta@f$, in rads) measured **counter-clockwise** when looking
    from the torus top, and an angle on the poloidal cross section (`w`, or
    `VMEC` @f$\theta@f$, also in rads). The covariant metric units are m^2, as
    required by the parent class. More info at the website
    [STELLOPT](https://princetonuniversity.github.io/STELLOPT/VMEC.html).
*/
class metric_vmec : public metric_covariant {
 public:
  typedef parser_vmec::narray_type narray_type;

  metric_vmec(
      const parser_vmec* parser, const interpolator1d_factory* ifactory);
  virtual ~metric_vmec() override {};

  virtual SM3 operator()(const IR3& position) const override;
  virtual dSM3 del(const IR3& position) const override;
  virtual IR3 transform2cylindrical(const IR3& position) const;

  const parser_vmec* parser() const { return parser_; };
 private:
  const size_t harmonics_;
  const narray_type m_, n_;
  const parser_vmec* parser_;
  std::vector<std::unique_ptr<interpolator1d>> r_mn_, z_mn_;

  typedef std::vector<std::complex<double>> cis_container_t;
  const cis_container_t& cached_cis(double theta, double zeta) const;
};

}  // end namespace gyronimo

#endif  // GYRONIMO_METRIC_VMEC
