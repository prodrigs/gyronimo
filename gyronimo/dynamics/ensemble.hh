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

// @ensemble.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_ENSEMBLE
#define GYRONIMO_ENSEMBLE

#include <ranges>
#include <execution>

namespace gyronimo {

//! Assembles a single `boost::odeint` system from a large `Gyron` collection.
/*!
    Takes a range of individual `Gyron` objects, each with a well defined
    operator `dfdt = gyron(f, t)`, and redefines a single collective operator
    `void ensemble(f, dfdt, t)` by running along the gyron collection and
    assembling the individual outputs. Notice that the class definition is
    independent of the particular container used by clent code to store the
    corresponding collection of individual states (i.e., `StateRange`).
*/
template<typename GyronRange>
class ensemble {
 public:
  typedef typename std::ranges::range_value_t<GyronRange> gyron_t;
  typedef typename gyron_t::state gyron_state_t;
  ensemble(const GyronRange& gyron_ensemble)
    : gyron_ensemble_(gyron_ensemble) {};
  template<typename StateRange>
  void operator()(
      const StateRange& f, StateRange& dfdt, double t) const {
    std::transform(std::execution::par,
        f.begin(), f.end(), gyron_ensemble_.begin(), dfdt.begin(),
            [&t](const gyron_state_t& s, const gyron_t& g) {return g(s, t);});
  };
  size_t size() const {return std::size(gyron_ensemble_);};
  auto begin() const {return std::begin(gyron_ensemble_);};
  auto end() const {return std::end(gyron_ensemble_);};
 private:
  const GyronRange& gyron_ensemble_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_ENSEMBLE
