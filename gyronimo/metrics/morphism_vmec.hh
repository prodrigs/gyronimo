// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Manuel Assunção.

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

// @morphism_vmec.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_VMEC
#define GYRONIMO_MORPHISM_VMEC

#include <gyronimo/interpolators/interpolator1d.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/parsers/parser_vmec.hh>

namespace gyronimo {

//! Morphism in `VMEC` curvilinear coordinates.
class morphism_vmec : public morphism {
 public:
  typedef std::valarray<double> narray_type;

  morphism_vmec(
      const parser_vmec* parser, const interpolator1d_factory* ifactory);
  virtual ~morphism_vmec() override;

  virtual IR3 operator()(const IR3& q) const override;
  virtual IR3 inverse(const IR3& x) const override;
  virtual IR3 translation(const IR3& q, const IR3& delta) const override;
  virtual dIR3 del(const IR3& q) const override;
  virtual ddIR3 ddel(const IR3& q) const override;

  virtual double jacobian(const IR3& q) const override;

  const parser_vmec* parser() const { return parser_; };
 private:
  const parser_vmec* parser_;
  narray_type xm_, xn_;

  interpolator1d** Rmnc_;
  interpolator1d** Zmns_;
  interpolator1d** gmnc_;

  std::pair<double, double> get_rz(const IR3& position) const;
  std::pair<double, double> reflection_past_axis(double s, double theta) const;
};

}  // end namespace gyronimo

#endif  // GYRONIMO_MORPHISM_VMEC
