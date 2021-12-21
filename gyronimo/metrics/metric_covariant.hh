// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @metric_covariant.hh

#ifndef GYRONIMO_METRIC_COVARIANT
#define GYRONIMO_METRIC_COVARIANT

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/core/SM3algebra.hh>

namespace gyronimo {

//! Abstract class implementing basic 3x3 metric-tensor functionality.
/*!
    Requires derived classes to implement the 6 independent components @f$
    g_{ij} @f$ of the covariant metric tensor for a given coordinate system in
    IR^3 [via the `operator()()` method] and their corresponding 18 partial
    derivatives [via the `del()` method]. The class implements general methods
    to compute the Jacobian (square root of its determinant) and its gradient,
    to transform contravariant vectors into covariant ones and conversely. These
    three methods are left as virtual to allow more efficient reimplementations
    in derived classes, if needed. The `position` argument stands for the three
    contravariant components @f$ q^i @f$ of the position vector @f$ \vec{q} @f$
    in each specific coordinate set. Regarding units, all methods must ensure
    that @f$g_{ij} q^i q^j@f$ returns values in SI (m^2).
*/
class metric_covariant {
 public:
  metric_covariant() {};
  virtual ~metric_covariant() {};

  virtual SM3 operator()(const IR3& r) const = 0;
  virtual dSM3 del(const IR3& r) const = 0;

  virtual double jacobian(const IR3& r) const;
  virtual IR3 del_jacobian(const IR3& r) const;
  virtual IR3 to_covariant(const IR3& B, const IR3& r) const;
  virtual IR3 to_contravariant(const IR3& B, const IR3& r) const;
  virtual SM3 inverse(const IR3& r) const;
  virtual dSM3 del_inverse(const IR3& r) const;
};

}

#endif // GYRONIMO_METRIC_COVARIANT