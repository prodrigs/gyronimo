// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues, Rogerio Jorge.

// @metric_stellna.hh

#ifndef GYRONIMO_METRIC_STELLNA
#define GYRONIMO_METRIC_STELLNA

#include <gyronimo/parsers/parser_stellna.hh>
#include <gyronimo/metrics/metric_covariant.hh>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo{

//! Covariant metric in stellarator near-axis coordinates.
/*!
    Supports indirect initialisation from a 'parser_stellna' object or direct
    initialisation from supplied arrays containing the solution of the
    Frenet-Serret equations for the magnetic axis sampled along the toroidal
    angle. The particular kind of 1d interpolator to use is set by 'ifactory'.

    The right-handed coordinate set @f$\{r, \theta, \phi\}@f$ is described in
    R.~Jorge et al., Nuc. Fusion **60**, 076021 (2020).
*/
class metric_stellna : public metric_covariant {
 public:
  metric_stellna(
      const parser_stellna *parser, const interpolator1d_factory* ifactory);
  metric_stellna(
      int field_periods, double eta_bar,
      const dblock& phi_grid, const dblock& sigma,
      const dblock& dldphi, const dblock& torsion, const dblock& curvature,
      const interpolator1d_factory* ifactory);
  virtual ~metric_stellna() override;

  virtual SM3 operator()(const IR3& position) const override;
  virtual dSM3 del(const IR3& position) const override;

  const interpolator1d* curvature() const {return curvature_;};
  const interpolator1d* torsion() const {return torsion_;};
  const interpolator1d* dldphi() const {return dldphi_;};
  const interpolator1d* sigma() const {return sigma_;};
  double reduce_phi(double phi) const;

 private:
  const double eta_bar_, phi_modulus_factor_;
  interpolator1d *sigma_, *curvature_, *torsion_, *dldphi_;
};

} // end namespace gyronimo

#endif // GYRONIMO_METRIC_STELLNA
