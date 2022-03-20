// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues.

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

// @field_line.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_FIELD_LINE
#define GYRONIMO_FIELD_LINE

#include <array>
#include <gyronimo/fields/IR3field.hh>

namespace gyronimo {

//! Field-line equations of motion parametrised by length over `Lref`.
class field_line {
 public:
  typedef std::array<double, 3> state;
  field_line(
      const IR3field* field, double Lref) : field_(field), Lref_(Lref) {};
  ~field_line() {};
  state operator()(const state& x, double t) const {
      return IR3(Lref_*field_->contravariant_versor(x, t));};
  const IR3field* field() const {return field_;};
  double Lref() const {return Lref_;};
 private:
  const IR3field* field_;
  const double Lref_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_FIELD_LINE
