// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues and Jorge Ferreira

// @metric_vmec.hh

#pragma once

#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/metrics/metric_covariant.hh>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo{

//! Covariant metric in `VMEC` curvilinear coordinates.
/*!
    Builds the metric from the information provided by a `parser_vmec` object.
*/
class metric_vmec : public metric_covariant {
 public:
  typedef std::valarray<double> narray_type;
  typedef std::vector<interpolator1d> spectralarray_type;

  metric_vmec(
    const parser_vmec *parser, const interpolator1d_factory *ifactory);
  virtual ~metric_vmec() override;
  virtual SM3 operator()(const IR3& position) const override;
  virtual dSM3 del(const IR3& position) const override;
  // virtual double jacobian(const IR3& position) const override;
  virtual IR3 transform2cylindrical(const IR3& position) const override;
  double jacobian_vmec(const IR3& position) const;
  const parser_vmec* parser() const {return parser_;};
  const double signgs() const {return signsgs_;};

 private:
  const parser_vmec* parser_;
  double b0_;
  int mnmax_, mnmax_nyq_, ns_, mpol_, ntor_, nfp_, signsgs_; 
  narray_type xm_, xn_, xm_nyq_, xn_nyq_; 
  // spectral interpolators
  interpolator1d **Rmnc_;
  interpolator1d **Zmns_;
  interpolator1d **gmnc_;
};

} // end namespace gyronimo
