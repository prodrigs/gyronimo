// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Jorge Ferreira and Paulo Rodrigues.

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

#include <gyronimo/interpolators/interpolator1d.hh>
#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism_vmec.hh>

namespace gyronimo {

//! Covariant metric corresponding to a `morphism_vmec` object.
/*
    This class inherits all functionality from its `metric_connected` parent,
    which is extracted from the `morphism_vmec` it is connected to.
*/
class metric_vmec : public metric_connected {
 public:
  metric_vmec(const morphism_vmec* morph);
  virtual ~metric_vmec() override {};
  const parser_vmec* my_parser() const { return parser_; };
  const morphism_vmec* my_morphism() const { return morphism_; };
 private:
  const morphism_vmec* morphism_;
  const parser_vmec* parser_;
};

}  // end namespace gyronimo

#endif  // GYRONIMO_METRIC_VMEC
