// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Jorge Ferreira and Paulo Rodrigues.

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

// @boris-xyz.cc, this file is part of ::gyronimo::

#include <gyronimo/dynamics/classical_boris.hh>
#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/metric_cartesian.hh>
#include <gyronimo/version.hh>

#include <iostream>
#include <cmath>

using namespace gyronimo;

void print_help() {
  std::cout << "boris-xyz, powered by ::gyronimo::v" << version_major << "."
            << version_minor << "." << version_patch
            << " (git-commit:" << git_commit_hash << ").\n";
  std::string help_message =
      "usage: boris-xyz x y z vx vy vz Bx By Bz qom time-step\n\n"
      "A trivial, command-line driven boris stepper over an homogeneous\n"
      "magnetic field in cartesian coordinates for pedagogical purposes.\n"
      "Takes the cartesian initial position and velocity, field components, \n"
      "charge over mass ratio, and time step from the command line (all in \n" 
      "SI units except qom, which is in proton charge/mass units) and prints \n"
      "the updated position and velocity to stdout.\n";
  std::cout << help_message;
  std::exit(0);
}

class homogeneneous_field : public IR3field {
 public:
  homogeneneous_field(const IR3& val, const metric_covariant* g)
      : IR3field(1, 1, g), val_(val) {};
  virtual ~homogeneneous_field() override {};
  virtual IR3 contravariant(const IR3& x, double time) const {return val_;};
  virtual IR3 covariant(const IR3& x, double time) const override {
    return val_;};
  virtual double magnitude(const IR3& x, double time) const override {
    return std::sqrt(inner_product(val_, val_));};
 private:
  IR3 val_;
};

std::array<double, 11> parse_args(char* argv[]) {
  return {
      std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3]),
      std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6]),
      std::atof(argv[7]), std::atof(argv[8]), std::atof(argv[9]),
      std::atof(argv[10]), std::atof(argv[11])};
}

int main(int argc, char* argv[]) {
  if (argc != 12) print_help();
  auto [x, y, z, vx, vy, vz, Bx, By, Bz, qom, time_step] = parse_args(argv);

  morphism_cartesian morphism;
  metric_cartesian g(&morphism);
  homogeneneous_field field({Bx, By, Bz}, &g);

  classical_boris boris(1, 1, qom, &field, nullptr);
  auto new_state = boris.do_step(
      boris.generate_state({x, y, z}, {vx, vy, vz}), 0.0, time_step);

  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  for (auto value : new_state) std::cout << value << " ";
  std:: cout << std::endl;

  return 0;
}
