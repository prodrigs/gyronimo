// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @io_streaming.hh

#ifndef GYRONIMO_IO_STREAMING
#define GYRONIMO_IO_STREAMING

#include <iostream>

namespace gyronimo {

//! Read helper for indexed objects.
template<class T>
void operator>>(std::istream& is,T& x) {
  for (size_t i = 0; i < x.size(); i++) is >> x[i];
}

//! Write helper for indexed objects.
template<class T>
void operator<<(std::ostream& os,T& x) {
  for (size_t i = 0; i < x.size(); i++) os << x[i];
}

} // end namespace gyronimo.

#endif // GYRONIMO_IO_STREAMING
