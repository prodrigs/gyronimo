// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @codata.hh

#ifndef GYRONIMO_CODATA
#define GYRONIMO_CODATA

namespace gyronimo {

//! Physical constants (SI units) by [CODATA](https://physics.nist.gov/cuu/Constants/index.html).
namespace codata {

constexpr double c = 2.99792458e8;
constexpr double e = 1.602176634e-19;
constexpr double mu0 = 1.25663706212e-6;
constexpr double m_alpha = 6.6446573357e-27;
constexpr double m_proton = 1.67262192369e-27;
constexpr double m_neutron = 1.67492749804e-27;
constexpr double m_electron = 9.1093837015e-31;

} // end namespace codata

} // end namespace gyronimo.

#endif // GYRONIMO_CODATA
