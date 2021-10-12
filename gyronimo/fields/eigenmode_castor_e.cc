// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @eigenmode_castor_e.cc

#include <gyronimo/fields/eigenmode_castor_e.hh>

namespace gyronimo{

IR3 eigenmode_castor_e::covariant(const IR3& position, double time) const {
  IR3 A = A_.covariant(position, time);
  return {std::real(-iw_*A[IR3::u]),
      std::real(-iw_*A[IR3::v]), std::real(-iw_*A[IR3::w])};
}
IR3 eigenmode_castor_e::contravariant(const IR3& position, double time) const {
  IR3 A = A_.contravariant(position, time);
  return {std::real(-iw_*A[IR3::u]),
      std::real(-iw_*A[IR3::v]), std::real(-iw_*A[IR3::w])};
}

} // end namespace gyronimo.
