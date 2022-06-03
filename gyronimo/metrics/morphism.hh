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

// @morphism.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM
#define GYRONIMO_MORPHISM

#include <gyronimo/core/IR3algebra.hh>

namespace gyronimo {

//! Abstract morphism from curvilinear (q) into cartesian (x) coordinates.
class morphism {
 public:
  virtual ~morphism() {};
  virtual IR3 operator()(const IR3& q) const = 0;
  virtual IR3 inverse(const IR3& x) const = 0;
  virtual dIR3 del(const IR3& q) const = 0;
  virtual dIR3 del_inverse(
      const IR3& q) const {return gyronimo::inverse(del(q));};
  virtual dIR3 tan_basis(const IR3& q) const {return del(q);};
  virtual dIR3 dual_basis(const IR3& q) const {return del_inverse(q);};
};

}

#endif // GYRONIMO_MORPHISM
