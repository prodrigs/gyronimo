// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues and Jorge Ferreira.

// @equilibrium_vmec.hh

#pragma once

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/metrics/metric_vmec.hh>
#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/interpolators/interpolator2d.hh>

namespace gyronimo{

//! Equilibrium magnetic field in 'VMEC' curvilinear coordinates.
/*!
    Following IR3Field rules, the magnetic field is normalised by a `m_factor`
    matching its on-axis value `B0()` in [T], which is located at `R0()` in [m].
    The coordinates are set by the `metric_vmec` object and the type of 1d
    interpolators is set by the specific `interpolator1d_factory` supplied.
    Contravariant components have dimensions of [m^{-1}]. Being an
    **equilibrium** field, `t_factor` is set to one.
    
    Only the minimal interface is implemented for the moment and further
    specialisations may enhance the object's performance.
*/
class equilibrium_vmec : public IR3field_c1{
 public:
  typedef std::valarray<double> narray_type;
  typedef std::vector<interpolator1d> spectralarray_type;
  equilibrium_vmec(
      const metric_vmec *g, const interpolator1d_factory *ifactory);
  virtual ~equilibrium_vmec() override;

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {return {0.0,0.0,0.0};};
  // double magnitude(const IR3& position, double time) const override;
  double magnitude_vmec(const IR3& position, double time) const;

  double R_0() const {return metric_->parser()->R_0();};
  double B_0() const {return metric_->parser()->B_0();};
  const metric_vmec* metric() const {return metric_;};
 private:
  const metric_vmec *metric_;
  narray_type xm_nyq_, xn_nyq_; 
  // spectral interpolators
  interpolator1d **bmnc_;
  interpolator1d **bsupumnc_;
  interpolator1d **bsupvmnc_;
};

}// end namespace gyronimo.