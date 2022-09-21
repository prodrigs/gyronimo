
#ifndef GYRONIMO_METRIC_CYLINDRICAL
#define GYRONIMO_METRIC_CYLINDRICAL

#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! Trivial covariant metric for cylindrical space.
class metric_cylindrical : public metric_covariant {
 public:
  metric_cylindrical() {};
  virtual ~metric_cylindrical() override {};

  virtual SM3 operator()(const IR3& q) const override {
    return {1.0, 0.0, 0.0, q[IR3::u]*q[IR3::u], 0.0, 1.0};
  };
  virtual dSM3 del(const IR3& q) const override {
    return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 2*q[IR3::u], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  };
  virtual double jacobian(const IR3& q) const override {
    return q[IR3::u];
  };
  virtual IR3 del_jacobian(const IR3& q) const override {
    return {std::cos(q[IR3::v]), std::sin(q[IR3::v]), 0.0};
  };
  virtual IR3 to_covariant(const IR3& B, const IR3& q) const override {
    return {B[IR3::u], B[IR3::v]*q[IR3::u]*q[IR3::u], B[IR3::w]};
  };
  virtual IR3 to_contravariant(const IR3& B, const IR3& q) const override {
    return {B[IR3::u], B[IR3::v]/(q[IR3::u]*q[IR3::u]), B[IR3::w]};
  };
}; // end class metric_cylindrical

} // end namespace gyronimo.

#endif // GYRONIMO_METRIC_CARTESIAN
