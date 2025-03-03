// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Jorge Ferreira and Paulo Rodrigues.

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

// @equilibrium_vmec_a.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_EQUILIBRIUM_VMEC_A
#define GYRONIMO_EQUILIBRIUM_VMEC_A

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/interpolators/interpolator2d.hh>
#include <gyronimo/metrics/metric_vmec.hh>

#include <complex>
#include <memory>

namespace gyronimo {

//! Equilibrium magnetic potential field in 'VMEC' curvilinear coordinates.
/*!
    This class represents the vector potential \( A \) of the equilibrium 
    magnetic field in VMEC coordinates. It follows the conventions of `IR3field`, 
    normalizing by setting `m_factor` to \( B_0 \) (in [T]) at \( R_0 \) (in [m]). 

    The vector potential is defined based on the `metric_vmec` object and 
    employs 1D interpolators provided by `interpolator1d_factory`. Contravariant 
    and covariant components are computed, along with their spatial derivatives. 

    Optional normalization of the field is controlled by the `normalised` flag.
*/
class equilibrium_a_field : public IR3field_c1 {
 public:
  using narray_type = parser_vmec::narray_type;
  equilibrium_a_field(
      const metric_vmec* g, const interpolator1d_factory* ifactory, bool normalised=true);
  virtual ~equilibrium_a_field() override {};

  virtual IR3 covariant(
      const IR3& position, double time) const override;
  virtual IR3 contravariant(
      const IR3& position, double time) const override;
  virtual dIR3 del_covariant(
      const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {return {0, 0, 0};};
  virtual IR3 partial_t_covariant(
      const IR3& position, double time) const override {return {0, 0, 0};};
  virtual double partial_t_magnitude(
      const IR3& position, double time) const override {return 0;};
  double magnitude_vmec(
      const IR3& position, double time) const;


  double R0() const { return parser_->R0(); };
  double B0() const { return parser_->B0(); };
  const metric_vmec* metric() const { return metric_; };
  const parser_vmec* my_parser() const { return parser_; };
  const morphism_vmec* my_morphism() const { return metric_->my_morphism(); };


 private:
  const metric_vmec* metric_;
  const parser_vmec* parser_;
  const size_t harmonics_;
  const narray_type m_, n_;
  std::vector<size_t> index_;
  std::vector<std::unique_ptr<interpolator1d>> btheta_mn_, bzeta_mn_, g_mn_;
  const bool normalised_;


  narray_type mc_, nc_;
  narray_type csupumnc, csupvmnc; 
  narray_type ctheta_integ, czeta_integ;
  std::vector<std::unique_ptr<interpolator1d>> ctheta_mn_, czeta_mn_;
  std::unique_ptr<interpolator1d> ctheta_00_, czeta_00_;
  void build_harmonics(narray_type& new_m, narray_type& new_n, 
                       const narray_type parser_m, const narray_type parser_n,
                       narray_type& c_theta_val, narray_type& c_zeta_val);

  void build_integrals(narray_type& ctheta_integ, narray_type& czeta_integ);



  void build_interpolator_array(
      std::vector<std::unique_ptr<interpolator1d>>& interpolator_array,
      const narray_type& samples_array, const interpolator1d_factory* ifactory, size_t size_harmonics);


  typedef std::vector<std::complex<double>> cis_container_t;
  const cis_container_t& cached_cis_new(double theta, double zeta) const;

};

}  // end namespace gyronimo.

#endif  // GYRONIMO_EQUILIBRIUM_A_FIELD
