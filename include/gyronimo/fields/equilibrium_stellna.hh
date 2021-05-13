// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues, Rogerio Jorge.

// @equilibrium_stellna.hh

#ifndef GYRONIMO_EQUILIBRIUM_STELLNA
#define GYRONIMO_EQUILIBRIUM_STELLNA

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/metrics/metric_stellna.hh>

namespace gyronimo{

//! Quasi-symmetric stellarator equilibrium field in near-axis coordinates.
/*!
    Following IR3Field rules, the magnetic field is normalised by a `m_factor`
    matching its on-axis value `axis_field()` in [T]. In turn, 'axis_length()'
    is in [m]. The coordinates are set by the `metric_stellna` object and the
    contravariant field components have dimensions of [m^{-1}]. Being an
    **equilibrium** field, `t_factor` is set to one.
*/
class equilibrium_stellna : public IR3field_c1{
 public:
  equilibrium_stellna(
      const metric_stellna *g,
      double axis_field, double axis_length, double axis_iota_minus_n);
  virtual ~equilibrium_stellna() override {};

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {return {0.0,0.0,0.0};};

  const metric_stellna* metric() const {return metric_;};
  double axis_field() const {return this->m_factor();};
  double axis_length() const {return axis_length_;};

 private:
  const metric_stellna *metric_;
  double axis_length_, length_factor_, iota_factor_;
};

}// end namespace gyronimo.

#endif // GYRONIMO_EQUILIBRIUM_STELLNA