// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @IR3fields_c1.hh

#ifndef GYRONIMO_IR3FIELD_C1
#define GYRONIMO_IR3FIELD_C1

#include <gyronimo/fields/IR3field.hh>

namespace gyronimo {

//! Time-dependent, continuously differentiable field in @f$\mathbb{R}^3@f$.
/*!
    Besides implementing the contravariant components of the normalised field at
    a position and time (inherited form base `IR3field`), derived classes must
    also implement their partial derivatives with respect to variables `IR3::u`,
    `IR3::v`, `IR3::w`, and time via the abstract methods del_contravariant()
    and partial_t_contravariant(). General purpose implementations of the curl
    operator (contravariant components returned), of the magnitude gradient
    (covariant components returned) and partial time derivative (scalar) are
    provided, as well as the partial derivatives of the covariant components
    [del_covariant() and partial_t_covariant()], which can be further optimized
    in derived classes if needed. Notice that time derivatives are taken with
    respect to the normalised time, whilst spacial derivatives are taken with
    respect to the arbitrary coordinates defined by the covariant_metric object.
    Conversion to SI values is achieved using the normalization constants
    `m_factor` and `t_factor`.
*/
class IR3field_c1 : public IR3field {
 public:
  IR3field_c1(double m_factor, double t_factor, const metric_covariant* g)
      : IR3field(m_factor, t_factor, g) {};
  virtual ~IR3field_c1() override {};

  virtual IR3 contravariant(const IR3& position, double time) const = 0;
  virtual dIR3 del_contravariant(const IR3& position, double time) const = 0;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const = 0;

  virtual IR3 del_magnitude(const IR3& position, double time) const;
  virtual double partial_t_magnitude(const IR3& position, double time) const;
  virtual dIR3 del_covariant( const IR3& position, double time) const;
  virtual IR3 partial_t_covariant(const IR3& position, double time) const;
  virtual IR3 curl(const IR3& position, double time) const;
};

} // end namespace gyronimo.

#endif // GYRONIMO_IR3FIELD_C1
