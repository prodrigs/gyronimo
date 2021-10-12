// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @equilibrium_helena.hh

#ifndef GYRONIMO_EQUILIBRIUM_HELENA
#define GYRONIMO_EQUILIBRIUM_HELENA

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/metrics/metric_helena.hh>
#include <gyronimo/parsers/parser_helena.hh>
#include <gyronimo/interpolators/interpolator2d.hh>

namespace gyronimo{

//! Tokamak equilibrium magnetic field in 'HELENA' curvilinear coordinates.
/*!
    Following IR3Field rules, the magnetic field is normalised by a `m_factor`
    matching its on-axis value `B0()` in [T], which is located at `R0()` in [m].
    The coordinates are set by the `metric_helena` object and the type of 2d
    interpolators is set by the specific `interpolator2d_factory` supplied.
    Contravariant components have dimensions of [m^{-1}]. Being an
    **equilibrium** field, `t_factor` is set to one.
    
    Only the minimal interface is implemented for the moment and further
    specialisations may enhance the object's performance.
*/
class equilibrium_helena : public IR3field_c1{
 public:
  equilibrium_helena(
      const metric_helena *g, const interpolator2d_factory *ifactory);
  virtual ~equilibrium_helena() override;

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {return {0.0,0.0,0.0};};

  double R0() const {return metric_->parser()->rmag();};
  double B0() const {return metric_->parser()->bmag();};

 private:
  const metric_helena *metric_;
  interpolator2d *Bchi_, *Bphi_;
};

}// end namespace gyronimo.


#endif // GYRONIMO_EQUILIBRIUM_HELENA
