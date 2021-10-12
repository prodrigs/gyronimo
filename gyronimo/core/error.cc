// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @error.cc

#include <cstdlib>
#include <iostream>
#include <gyronimo/core/error.hh>

namespace gyronimo {

//! Prints a warning message to `stderr`.
void warning(const char* message) {
  std::cerr << "gyronimo::warning: " << message << std::endl;
}

//! Prints an error message to `stderr` and exits execution with `exit_code`.
void error(
    const char* caller_name, const char* file_name,
    const int line_num, const char* message, const int exit_code) {
  std::cerr << "gyronimo::" << caller_name << ": " << message << std::endl;
  std::cerr << file_name << ", line " << line_num << ".\n";
  std::cerr << "error code: " << exit_code << ".\n";
  std::exit(exit_code);
}

} // end namespace gyronimo.
