// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @metric_stellna.hh

#ifndef GYRONIMO_METRIC_STELLNA
#define GYRONIMO_METRIC_STELLNA

#include <gyronimo/parsers/parser_stellna.hh>
#include <gyronimo/metrics/metric_covariant.hh>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo{

//! Covariant metric in stellarator near-axis coordinates.
class metric_stellna : public metric_covariant {
 public:
  metric_stellna(
      const parser_stellna *parser, const interpolator1d_factory *ifactory);
  virtual ~metric_stellna() override;

  virtual SM3 operator()(const IR3& position) const override;
  virtual dSM3 del(const IR3& position) const override;

  const parser_stellna* parser() const {return parser_;};

 private:
  const parser_stellna* parser_;
  const double phi_modulus_factor_;
  interpolator1d *sigma_, *curvature_, *torsion_, *dldphi_;
};

} // end namespace gyronimo

#endif // GYRONIMO_METRIC_STELLNA
