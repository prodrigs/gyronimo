// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Jorge Ferreira and Paulo Rodrigues.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @metric_vmec.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_VMEC
#define GYRONIMO_METRIC_VMEC

#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/metrics/morphism_vmec.hh>
#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/interpolators/interpolator1d.hh>

namespace gyronimo{

//! Covariant metric in `VMEC` curvilinear coordinates.
class metric_vmec : public metric_connected {
 public:
  typedef std::valarray<double> narray_type;
  typedef std::vector<interpolator1d> spectralarray_type;

  metric_vmec(const morphism_vmec* morph, 
    const interpolator1d_factory *ifactory);
  virtual ~metric_vmec() override;
  virtual SM3 operator()(const IR3& position) const override;
  virtual dSM3 del(const IR3& position) const override;

  const parser_vmec* parser() const {return parser_;};
  const morphism_vmec* my_morphism() const {
    return static_cast<const morphism_vmec*>(metric_connected::my_morphism());
  };

 private:
  const parser_vmec* parser_;
  const interpolator1d_factory *ifactory_;
  int mnmax_, mnmax_nyq_, ns_, mpol_, ntor_, nfp_, signsgs_; 
  narray_type xm_, xn_, xm_nyq_, xn_nyq_; 
  interpolator1d **Rmnc_;
  interpolator1d **Zmns_;
  interpolator1d **gmnc_;
};

} // end namespace gyronimo

#endif // GYRONIMO_METRIC_VMEC
