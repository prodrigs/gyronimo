// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @error.hh

#ifndef GYRONIMO_ERROR
#define GYRONIMO_ERROR

namespace gyronimo {

void warning(const char* message);
void error(
    const char* caller_name, const char* file_name,
    const int line_num, const char* message, const int exit_code);

} // end namespace gyronimo.

#endif // GYRONIMO_ERROR.
