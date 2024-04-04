// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Manuel Assunção and Paulo Rodrigues.

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

// @morphism_vmec.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_VMEC
#define GYRONIMO_MORPHISM_VMEC

#include <gyronimo/interpolators/interpolator1d.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/parsers/parser_vmec.hh>

#include <complex>
#include <memory>
#include <numbers>

namespace gyronimo {

//! Morphism in `VMEC` curvilinear coordinates.
/*!
    Builds the morphism data structure from the information provided by a
    `parser_vmec` object. The right-handed coordinates are the poloidal flux
    normalised to its boundary value (`u`, or `VMEC` @f$\chi@f$), the toroidal
    angle (`v`, or `VMEC` @f$\zeta@f$, in rads) measured **counter-clockwise**
    when looking from the torus top, and an angle on the poloidal cross section
    (`w`, or `VMEC` @f$\theta@f$, also in rads). More info at the website
    [STELLOPT](https://princetonuniversity.github.io/STELLOPT/VMEC.html).
*/
class morphism_vmec : public morphism {
 public:
  using narray_type = parser_vmec::narray_type;
  morphism_vmec(
      const parser_vmec* parser, const interpolator1d_factory* ifactory);
  virtual ~morphism_vmec() override {};
  virtual IR3 operator()(const IR3& q) const override;
  virtual IR3 inverse(const IR3& x) const override;
  virtual double jacobian(const IR3& q) const override;
  virtual dIR3 del(const IR3& q) const override;
  virtual ddIR3 ddel(const IR3& q) const override;
  virtual IR3 translation(const IR3& q, const IR3& delta) const override;

  const parser_vmec* my_parser() const { return parser_; };
  std::pair<double, double> get_rz(const IR3& q) const;
 private:
  const parser_vmec* parser_;
  const size_t harmonics_;
  const narray_type m_, n_;
  std::vector<size_t> index_;
  std::vector<std::unique_ptr<interpolator1d>> r_mn_, z_mn_;

  using cis_container_t = std::vector<std::complex<double>>;
  const cis_container_t& cached_cis(double theta, double zeta) const;
  IR3 inverse(const IR3& x, const std::pair<double, double>& guess) const;
  std::pair<double, double> reflection_past_axis(
      double flux, double theta) const;
  void build_interpolator_array(
      std::vector<std::unique_ptr<interpolator1d>>& interpolator_array,
      const narray_type& samples_array, const interpolator1d_factory* ifactory);
  struct aux_rz_t {
    double r, z;
  };
  struct aux_rz_del_t {
    double r, z, drdu, drdw, dzdu, dzdw;
  };
  struct aux_del_t {
    double r, drdu, drdv, drdw, dzdu, dzdv, dzdw;
  };
  struct aux_ddel_t {
    double r, drdu, drdv, drdw, dzdu, dzdv, dzdw;
    double d2rdudu, d2rdudv, d2rdudw, d2rdvdv, d2rdvdw, d2rdwdw;
    double d2zdudu, d2zdudv, d2zdudw, d2zdvdv, d2zdvdw, d2zdwdw;
  };
  friend aux_rz_t operator+(const aux_rz_t& x, const aux_rz_t& y);
  friend aux_rz_del_t operator+(const aux_rz_del_t& x, const aux_rz_del_t& y);
  friend aux_del_t operator+(const aux_del_t& x, const aux_del_t& y);
  friend aux_ddel_t operator+(const aux_ddel_t& x, const aux_ddel_t& y);
};

inline IR3 morphism_vmec::operator()(const IR3& q) const {
  double zeta = q[IR3::v];
  auto [r, z] = this->get_rz(q);
  return {r * std::cos(zeta), r * std::sin(zeta), z};
}

inline IR3 morphism_vmec::translation(const IR3& q, const IR3& delta) const {
  return this->inverse((*this)(q) + delta, {q[IR3::u], q[IR3::w]});
}

inline std::pair<double, double> morphism_vmec::reflection_past_axis(
    double flux, double theta) const {
  return flux < 0 ?
      std::pair<double, double> {-flux, theta + std::numbers::pi} :
      std::pair<double, double> {flux, theta};
}

inline morphism_vmec::aux_rz_t operator+(
    const morphism_vmec::aux_rz_t& x, const morphism_vmec::aux_rz_t& y) {
  return {x.r + y.r, x.z + y.z};
}

inline morphism_vmec::aux_rz_del_t operator+(
    const morphism_vmec::aux_rz_del_t& x,
    const morphism_vmec::aux_rz_del_t& y) {
  return {
      x.r + y.r, x.z + y.z, x.drdu + y.drdu, x.drdw + y.drdw, x.dzdu + y.dzdu,
      x.dzdw + y.dzdw};
}

inline morphism_vmec::aux_del_t operator+(
    const morphism_vmec::aux_del_t& x, const morphism_vmec::aux_del_t& y) {
  return {
      x.r + y.r, x.drdu + y.drdu, x.drdv + y.drdv, x.drdw + y.drdw,
      x.dzdu + y.dzdu, x.dzdv + y.dzdv, x.dzdw + y.dzdw};
}

inline morphism_vmec::aux_ddel_t operator+(
    const morphism_vmec::aux_ddel_t& x, const morphism_vmec::aux_ddel_t& y) {
  return {
      x.r + y.r, x.drdu + y.drdu, x.drdv + y.drdv, x.drdw + y.drdw,
      x.dzdu + y.dzdu, x.dzdv + y.dzdv, x.dzdw + y.dzdw, x.d2rdudu + y.d2rdudu,
      x.d2rdudv + y.d2rdudv, x.d2rdudw + y.d2rdudw, x.d2rdvdv + y.d2rdvdv,
      x.d2rdvdw + y.d2rdvdw, x.d2rdwdw + y.d2rdwdw, x.d2zdudu + y.d2zdudu,
      x.d2zdudv + y.d2zdudv, x.d2zdudw + y.d2zdudw, x.d2zdvdv + y.d2zdvdv,
      x.d2zdvdw + y.d2zdvdw, x.d2zdwdw + y.d2zdwdw};
}

}  // end namespace gyronimo

#endif  // GYRONIMO_MORPHISM_VMEC
