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

// @vmectrace.cc, this file is part of ::gyronimo::

// Command-line tool to print guiding-centre orbits in `VMEC` equilibria.
// External dependencies:
// - [argh](https://github.com/adishavit/argh), a minimalist argument handler.
// - [GSL](https://www.gnu.org/software/gsl), the GNU Scientific Library.
// - [boost](https://www.boost.org), the boost library.
// - [netcdf-c++4] (https://github.com/Unidata/netcdf-cxx4.git).

#include <gyronimo/core/codata.hh>
#include <gyronimo/dynamics/guiding_centre.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>
#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/version.hh>

#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <argh.h>
#include <cmath>
#include <iostream>

using namespace gyronimo;

void print_help() {
  std::cout << "vmectrace, powered by ::gyronimo::v" << version_major
            << "." << version_minor << "." << version_patch << ".\n";
  std::string help_message =
      "usage: vmectrace [options] vmec_netcdf_file\n"
      "reads a vmec output file, prints guiding-centre orbit to stdout.\n"
      "options:\n"
      "  -lref= Reference length (in si, default 1).\n"
      "  -vref= Reference velocity (in si, default 1).\n"
      "  -flux= Initial toroidal flux (vmec, default 0.5).\n"
      "  -mass= Particle mass (in m_proton, default 1).\n"
      "  -charge=\n"
      "         Particle charge (in q_proton, default 1).\n"
      "  -zeta=, -theta=\n"
      "         Initial zeta and theta (vmec angles in rad, default 0).\n"
      "  -energy=, -lambda=\n"
      "         Energy (eV) and lambda signed as v_parallel (default 1).\n"
      "  -tfinal=, -samples=\n"
      "         Time limit (lref/vref, default 1) and samples (default 512).\n"
      "  Note: lambda=magnetic_moment_si*B_axis_si/energy_si.\n";
  std::cout << help_message;
  std::exit(0);
}

class orbit_observer {
 public:
  orbit_observer(const equilibrium_vmec* e, const guiding_centre* g)
      : eq_pointer_(e), gc_pointer_(g) {};
  void operator()(const guiding_centre::state& s, double t) {
    IR3 q = gc_pointer_->get_position(s);
    auto [R, z] = eq_pointer_->my_morphism()->get_rz(q);
    double phi = q[IR3::v], x = R * std::cos(phi), y = R * std::sin(phi);
    std::cout << t << " " << q[IR3::u] << " " << q[IR3::v] << " " << q[IR3::w]
              << " " << gc_pointer_->energy_perpendicular(s, t) << " "
              << gc_pointer_->energy_parallel(s) << " " << x << " " << y << " "
              << z << '\n';
  };
 private:
  const equilibrium_vmec* eq_pointer_;
  const guiding_centre* gc_pointer_;
};

int main(int argc, char* argv[]) {
  auto command_line = argh::parser(argv);
  if (command_line[{"h", "help"}]) print_help();
  if (!command_line(1)) {  // the 1st non-option argument is the mapping file.
    std::cout << "vmectrace: no vmec equilibrium file provided; -h for help.\n";
    std::exit(1);
  }
  cubic_gsl_factory ifactory;
  parser_vmec parser(command_line[1]);
  gyronimo::morphism_vmec morph(&parser, &ifactory);
  gyronimo::metric_vmec g(&morph);
  equilibrium_vmec veq(&g, &ifactory);

  double flux, zeta, mass, lref, vref, theta, tfinal, charge, energy, lambda;
  command_line("flux", 0.5) >> flux;
  command_line("zeta", 0.0) >> zeta;
  command_line("mass", 1.0) >> mass;
  command_line("lref", 1.0) >> lref;
  command_line("vref", 1.0) >> vref;
  command_line("theta", 0.0) >> theta;
  command_line("tfinal", 1.0) >> tfinal;
  command_line("charge", 1.0) >> charge;
  command_line("energy", 1.0) >> energy;
  command_line("lambda", 1.0) >> lambda;
  double vpp_sign = std::copysign(1.0, lambda);  // lambda carries vpp sign.
  lambda = std::abs(lambda);  // once vpp sign is stored, lambda turns unsigned.

  double energy_ref = 0.5 * codata::m_proton * mass * vref * vref;
  double energy_si = energy * codata::e;
  guiding_centre gc(
      lref, vref, charge / mass, lambda * energy_si / energy_ref, &veq, nullptr);
  guiding_centre::state initial_state = gc.generate_state(
      {flux, zeta, theta}, energy_si / energy_ref,
      (vpp_sign > 0 ? guiding_centre::plus : guiding_centre::minus), 0);

  std::cout << "# vmectrace, powered by ::gyronimo::v" << version_major << "."
            << version_minor << " " << version_patch << ".\n";
  std::cout << "# args: ";
  for (int i = 1; i < argc; i++) std::cout << argv[i] << " ";
  std::cout << std::endl
            << "# E_ref: " << energy_ref << " [J]"
            << " B_axis: " << veq.m_factor() << " [T]"
            << " mu_tilde: " << gc.mu_tilde() << '\n';
  std::cout << "# vars: t flux zeta theta E_perp/E_ref E_parallel/E_ref x y "
               "z\n";

  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  size_t nsamples;
  command_line("samples", 512) >> nsamples;
  boost::numeric::odeint::runge_kutta4<guiding_centre::state>
      integration_algorithm;
  boost::numeric::odeint::integrate_const(
      integration_algorithm, odeint_adapter(&gc), initial_state, 0.0, tfinal,
      tfinal / nsamples, orbit_observer(&veq, &gc));

  return 0;
}
