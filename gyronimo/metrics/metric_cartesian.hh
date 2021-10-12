// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @metric_cartesian.hh

#ifndef GYRONIMO_METRIC_CARTESIAN
#define GYRONIMO_METRIC_CARTESIAN

#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! Trivial covariant metric for cartesian space.
class metric_cartesian : public metric_covariant {
 public:
  metric_cartesian() {};
  virtual ~metric_cartesian() override {};

  virtual SM3 operator()(const IR3& r) const override {
    return {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  };
  virtual dSM3 del(const IR3& r) const override {
    return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  };
  virtual double jacobian(const IR3& r) const override {
    return 1.0;
  };
  virtual IR3 del_jacobian(const IR3& r) const override {
    return {0.0, 0.0, 0.0};
  };
  virtual IR3 to_covariant(const IR3& B, const IR3& r) const override {
    return B;
  };
  virtual IR3 to_contravariant(const IR3& B, const IR3& r) const override {
    return B;
  };
};

} // end namespace gyronimo.

#endif // GYRONIMO_METRIC_CARTESIAN
