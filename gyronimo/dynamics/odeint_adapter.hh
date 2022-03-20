// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @odeint_adapter.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_ODEINT_ADAPTER
#define GYRONIMO_ODEINT_ADAPTER

namespace gyronimo {

//! Adapts a class F to work as an `ODEint` dynamical system.
//  @todo introduce concepts  to check for `F::operator()` and `F::state`.
template<class F>
class odeint_adapter {
 public:
  odeint_adapter(const F* g) : p_(g) {};
  void operator()(
      const F::state& x, F::state& dxdt, double t) const {dxdt = (*p_)(x, t);};
 private:
  const F* p_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_ODEINT_ADAPTER

