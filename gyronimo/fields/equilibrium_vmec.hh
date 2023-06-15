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

// @equilibrium_vmec.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_EQUILIBRIUM_VMEC
#define GYRONIMO_EQUILIBRIUM_VMEC

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/interpolators/interpolator2d.hh>
#include <gyronimo/metrics/metric_vmec.hh>
#include <gyronimo/parsers/parser_vmec.hh>

namespace gyronimo {

//! Equilibrium magnetic field in 'VMEC' curvilinear coordinates.
/*!
    Following `IR3field` rules, the magnetic field is normalised by setting
    `m_factor` to its value `B0()` (in [T]) on axis, which is located at `R0()`
    (in [m]). The coordinates are defined by the `metric_vmec` object and the
    type of 1d interpolators is set by the specific `interpolator1d_factory`
    supplied. Contravariant components have dimensions of [m^{-1}]. Being an
    **equilibrium** field, `t_factor` is set to one. Only the minimal interface
    is implemented here, all other functionality is inherited from the parent
    classes.
*/
class equilibrium_vmec : public IR3field_c1 {
 public:
  typedef parser_vmec::narray_type narray_type;
  equilibrium_vmec(
      const metric_vmec* g, const interpolator1d_factory* ifactory);
  virtual ~equilibrium_vmec() override {};

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {
    return {0.0, 0.0, 0.0};
  };

  double R0() const { return metric_->parser()->R0(); };
  double B0() const { return metric_->parser()->B0(); };
  const metric_vmec* metric() const { return metric_; };
 private:
  const size_t harmonics_;
  const narray_type m_, n_;
  const metric_vmec* metric_;
  std::vector<std::unique_ptr<interpolator1d>> btheta_mn_, bzeta_mn_;

  typedef std::vector<std::complex<double>> cis_container_t;
  const cis_container_t& cached_cis(double theta, double zeta) const;
};

}  // end namespace gyronimo.

#endif  // GYRONIMO_EQUILIBRIUM_VMEC
