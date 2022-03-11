// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @metric_spherical.hh

#ifndef GYRONIMO_METRIC_SPHERICAL
#define GYRONIMO_METRIC_SPHERICAL

#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! Covariant metric for spherical coordinates.
/*!
    The three contravariant coordinates are the distance to the origin
    normalized to `radius_norm` (`u`, with `radius_norm` in SI), the polar angle
    measured from the z-axis (co-latitude `v`, in rads), and the azimuthal angle
    (`w`, also in rads) measured clockwise when looking from the origin along
    the z-axis. Some inherited methods are overriden for efficiency.
*/
class metric_spherical : public metric_covariant {
 public:
  metric_spherical(double radius_norm);
  virtual ~metric_spherical() override {};

  virtual SM3 operator()(const IR3& r) const override ;
  virtual dSM3 del(const IR3& r) const override;

  virtual double jacobian(const IR3& r) const override;
  virtual IR3 del_jacobian(const IR3& r) const override;
  virtual IR3 to_covariant(const IR3& B, const IR3& r) const override ;
  virtual IR3 to_contravariant(const IR3& B, const IR3& r) const override;
 private:
  const double radius_norm_;
  const double radius_norm_squared_, radius_norm_cube_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_METRIC_SPHERICAL
