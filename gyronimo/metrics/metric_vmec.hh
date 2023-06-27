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

#include <gyronimo/interpolators/interpolator1d.hh>
#include <gyronimo/metrics/metric_covariant.hh>
#include <gyronimo/parsers/parser_vmec.hh>

#include <complex>
#include <memory>

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
  std::vector<size_t> index_;
  std::vector<std::unique_ptr<interpolator1d>> r_mn_, z_mn_;

  void build_interpolator_array(
      std::vector<std::unique_ptr<interpolator1d>>& interpolator_array,
      const narray_type& samples_array, const interpolator1d_factory* ifactory);

  typedef std::vector<std::complex<double>> cis_container_t;
  const cis_container_t& cached_cis(double theta, double zeta) const;

  struct auxiliar1_t {
    double r, drdu, drdv, drdw, dzdu, dzdv, dzdw;
  };
  friend auxiliar1_t operator+(const auxiliar1_t& x, const auxiliar1_t& y);
  struct auxiliar2_t {
    double r, z, drdu, drdv, drdw, dzdu, dzdv, dzdw;
    double d2rdudu, d2rdudv, d2rdudw, d2rdvdv, d2rdvdw, d2rdwdw;
    double d2zdudu, d2zdudv, d2zdudw, d2zdvdv, d2zdvdw, d2zdwdw;
  };
  friend auxiliar2_t operator+(const auxiliar2_t& x, const auxiliar2_t& y);
};

inline metric_vmec::auxiliar1_t operator+(
    const metric_vmec::auxiliar1_t& x, const metric_vmec::auxiliar1_t& y) {
  return {x.r + y.r, x.drdu + y.drdu, x.drdv + y.drdv, x.drdw + y.drdw,
          x.dzdu + y.dzdu, x.dzdv + y.dzdv, x.dzdw + y.dzdw};
}

inline metric_vmec::auxiliar2_t operator+(
    const metric_vmec::auxiliar2_t& x, const metric_vmec::auxiliar2_t& y) {
  return {x.r + y.r, x.z + y.z,
      x.drdu + y.drdu, x.drdv + y.drdv, x.drdw + y.drdw,
      x.dzdu + y.dzdu, x.dzdv + y.dzdv, x.dzdw + y.dzdw,
      x.d2rdudu + y.d2rdudu, x.d2rdudv + y.d2rdudv, x.d2rdudw + y.d2rdudw,
      x.d2rdvdv + y.d2rdvdv, x.d2rdvdw + y.d2rdvdw, x.d2rdwdw + y.d2rdwdw,
      x.d2zdudu + y.d2zdudu, x.d2zdudv + y.d2zdudv, x.d2zdudw + y.d2zdudw,
      x.d2zdvdv + y.d2zdvdv, x.d2zdvdw + y.d2zdvdw, x.d2zdwdw + y.d2zdwdw};
}

}  // end namespace gyronimo

#endif  // GYRONIMO_METRIC_VMEC
