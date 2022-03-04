// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Rogerio Jorge, Paulo Rodrigues.

// @equilibrium_stellnaqs.hh

#ifndef GYRONIMO_EQUILIBRIUM_STELLNAQS
#define GYRONIMO_EQUILIBRIUM_STELLNAQS

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/metrics/metric_stellnaqs.hh>

namespace gyronimo{

//! Quasi-symmetric stellarator equilibrium field via near-axis expansion.
/*!
    Following `IR3Field` rules, the magnetic field is normalised by the its
    on-axis value `axis_field` in [T] whereas `axis_length` is in [m]. The
    coordinates are set by the `metric_stellnaqs` object and the contravariant
    field components have dimensions of [m^{-1}]. Being an **equilibrium**
    field, `t_factor` is set to one.
*/
class equilibrium_stellnaqs : public IR3field_c1{
 public:
  equilibrium_stellnaqs(
      const metric_stellnaqs *g,
      double axis_field, double axis_length, double axis_iota);
  virtual ~equilibrium_stellnaqs() override {};

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {return {0.0,0.0,0.0};};

  const metric_stellnaqs* metric() const {return metric_;};
  double axis_field() const {return this->m_factor();};
  double axis_length() const {return axis_length_;};

 private:
  const metric_stellnaqs *metric_;
  double axis_length_, length_factor_, iota_factor_;
};

}// end namespace gyronimo.

#endif // GYRONIMO_EQUILIBRIUM_STELLNAQS
