// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @akima_boost.cc

#include <gyronimo/core/error.hh>
#include <gyronimo/interpolators/akima_boost.hh>

namespace gyronimo {

double akima_boost::derivative2(double x) const {
  error(__func__, __FILE__, __LINE__, "prime2 not available", 1);
  return 0.0;
}

interpolator1d* akima_boost_factory::interpolate_data(
    const dblock& abcissas, const dblock& ordinates) const {
  if (abcissas.size() != ordinates.size())
      error(__func__, __FILE__, __LINE__, "size mismatch", 1);
  switch(policy_) {
    case native:
      return new akima_boost(abcissas, ordinates);
      break;
    case periodic: {
      size_t n = std::size(ordinates);
      double double_step =
          (abcissas[1] - abcissas[0]) + (abcissas[n - 1] - abcissas[n - 2]);
      double left_right_prime =
          (ordinates[1] - ordinates[n - 2])/double_step;
      return new akima_boost(
          abcissas, ordinates, left_right_prime, left_right_prime);
      break;
    }
    default: error(__func__, __FILE__, __LINE__, "unknown policy", 1);
  }
  return nullptr;
}

} // end namespace gyronimo.
