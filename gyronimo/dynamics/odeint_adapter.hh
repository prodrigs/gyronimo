// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @odeint_adapter.hh

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

