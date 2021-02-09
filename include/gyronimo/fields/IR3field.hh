// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @IR3field.cc

#ifndef GYRONIMO_IR3FIELD
#define GYRONIMO_IR3FIELD

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! Base class for *adimensional* time-dependent fields in @f$\mathbb{R}^3@f$.
/*!
    Normalisation rules:
    1. **only** adimensional fields can be represented;
    2. Physical units can be restored by the magnitude factor `m_factor`;
    3. time is normalised to `t_factor` and also adimensional;
    4. `m_factor` and `t_factor` are provided to the constructor in SI units.
    
    Derived classes must implement the virtual method `contravariant(IR3&,
    double)`, returning the contravariant components of the adimensional field
    at a given position and time. The position coordinates are defined by the
    specific `metric_covariant` object supplied to the constructor, which also
    sets the physical dimensions of the contravariant components.

    The remaining interface [i.e., `covariant(...)`, `magnitude(...)`,
    `covariant_versor(...)`, and `contravariant_versor(...)`] is provided in
    general terms, but is left virtual to allow optimized implementations in
    derived classes.
*/
class IR3field {
 public:
  IR3field(double m_factor, double t_factor, const metric_covariant* g);
  IR3field(const IR3field& field)
      : m_factor_(field.m_factor_), t_factor_(field.t_factor_), 
        covariant_metric_(field.covariant_metric_) {};
  virtual ~IR3field() {};

  virtual IR3 contravariant(const IR3& position, double time) const = 0;

  virtual IR3 covariant(const IR3& position, double time) const;
  virtual double magnitude(const IR3& position, double time) const;
  virtual IR3 covariant_versor(const IR3& position, double time) const;
  virtual IR3 contravariant_versor(const IR3& position, double time) const;

  double m_factor() const {return m_factor_;};
  double t_factor() const {return t_factor_;};
  const metric_covariant* metric() const {return covariant_metric_;};

 private:
  const double m_factor_, t_factor_;
  const metric_covariant *covariant_metric_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_IR3FIELD
