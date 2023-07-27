// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2023 Paulo Rodrigues.

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

// @predicated_gyron.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_PREDICATED_GYRON
#define GYRONIMO_PREDICATED_GYRON

#include <functional>

namespace gyronimo {

//! Overrides a Gyron's `operator()` with `value` if `predicate` yields true.
template<typename Gyron>
class predicated_gyron {
 public:
  typedef Gyron::state state_t;
  typedef std::function<bool(const state_t)> predicate_t;
  template<typename... Args>
  predicated_gyron(predicate_t&& predicate, state_t& value, Args... gyron_args)
      : gyron_(gyron_args...), predicate_val_(value), predicate_(predicate) {};
  state_t operator()(const state_t& s, double t) const {
    return predicate_(s) ? gyron_(s, t) : predicate_val_;
  };
  const Gyron& get() const { return gyron_; };
 private:
  const Gyron gyron_;
  const state_t predicate_val_;
  const predicate_t predicate_;
};

}  // end namespace gyronimo.

#endif  // GYRONIMO_PREDICATED_GYRON
