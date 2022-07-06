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

// @vmectrace.cc, this file is part of ::gyronimo::

// Command-line tool to print guiding-centre orbits in `VMEC` equilibria.
// External dependencies:
// - [argh](https://github.com/adishavit/argh), a minimalist argument handler.
// - [GSL](https://www.gnu.org/software/gsl), the GNU Scientific Library.
// - [boost](https://www.gnu.org/software/gsl), the boost library.

#include <cmath>
#include <iostream>

#include <argh.h>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>

#include <gyronimo/version.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>
#include <gyronimo/dynamics/guiding_centre.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

void print_help() {
  std::cout << "vmectrace, powered by ::gyronimo::"
      << gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
  std::cout <<
      "usage: vmectrace [options] vmap\n";
  std::cout <<
      "reads an VMEC vmap and prints the required orbit to stdout.\n";
  std::cout << "options:\n";
  std::cout << "  -lref=val      Reference length (in si, default 1).\n";
  std::cout << "  -vref=val      Reference velocity (in si, default 1).\n";
  std::cout << "  -mass=val      Particle mass (in m_proton, default 1).\n";
  std::cout << "  -rho=val       Initial radius (vmec, default 0.5).\n";
  std::cout << "  -zeta=val      Initial zeta (vmec in rad, default 0).\n";
  std::cout << "  -theta=val     Initial theta (vmec in rad, default 0).\n";
  std::cout << "  -energy=val    Energy value (in eV, default 1).\n";
  std::cout << "  -lambda=val    Lambda value, signed as v_par (default 1).\n";
  std::cout << "  -tfinal=val    Time limit (in lref/vref, default 1).\n";
  std::cout << "  -charge=val    Particle charge (in q_proton, default 1).\n";
  std::cout << "  -nsamples=val  Number of orbit samples (default 512).\n";
  std::cout << "Note: lambda=magnetic_moment_si*B_axis_si/energy_si.\n";
  std::exit(0);
}

// ODEInt observer object to print diagnostics at each time step.
using namespace gyronimo;
class orbit_observer {
 public:
  orbit_observer(
      const IR3field_c1* e, const guiding_centre* g)
    : eq_pointer_(e), gc_pointer_(g) {};
  void operator()(const guiding_centre::state& s, double t) {
    IR3 q = gc_pointer_->get_position(s);
    IR3 c = eq_pointer_->metric()->transform2cylindrical(q);
    double R = c[IR3::u], phi = c[IR3::v], z = c[IR3::w];
    double x = R*std::cos(phi), y = R*std::sin(phi);
    std::cout << t << " "
        << q[IR3::u] << " "
        << q[IR3::v] << " "
        << q[IR3::w] << " "
        << gc_pointer_->energy_perpendicular(s, t) << " "
        << gc_pointer_->energy_parallel(s) << " " 
        << x << " " << y << " " << z << "\n";
  };
 private:
  const IR3field_c1* eq_pointer_;
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
  metric_vmec g(&parser, &ifactory);
  equilibrium_vmec veq(&g, &ifactory);

// Reads parameters from the command line:
  double rho; command_line("rho", 0.5) >> rho;
  double zeta; command_line("zeta", 0.0) >> zeta;
  double mass; command_line("mass", 1.0) >> mass;
  double lref; command_line("lref", 1.0) >> lref;
  double vref; command_line("vref", 1.0) >> vref;
  double theta; command_line("theta", 0.0) >> theta;
  double tfinal; command_line("tfinal", 1.0) >> tfinal;
  double charge; command_line("charge", 1.0) >> charge;
  double energy; command_line("energy", 1.0) >> energy;
  double lambda; command_line("lambda", 1.0) >> lambda;
  double vpp_sign = std::copysign(1.0, lambda);  // lambda carries vpp sign.
  lambda = std::abs(lambda);  // once vpp sign is stored, lambda turns unsigned.

  double energy_ref = 0.5*codata::m_proton*mass*vref*vref;
  double energy_si = energy*codata::e;
  guiding_centre gc(lref, vref, charge/mass, lambda*energy_si/energy_ref, &veq);
  guiding_centre::state initial_state = gc.generate_state(
      {rho, zeta, theta}, energy_si/energy_ref,
          (vpp_sign > 0 ?  guiding_centre::plus : guiding_centre::minus));

// Prints output header:
  std::cout << "# vmectrace, powered by ::gyronimo::"
      << version_major << "." << version_minor << ".\n";
  std::cout << "# args: ";
  for(int i = 1; i < argc; i++) std::cout << argv[i] << " ";
  std::cout << std::endl
      << "# E_ref: " << energy_ref << " [J]"
          << " B_axis: " << veq.m_factor() << " [T]"
              << " mu_tilde: " << gc.mu_tilde() << "\n";
  std::cout << "# vars: t rho zeta theta E_perp/E_ref E_parallel/E_ref x y z\n";

// integrates for t in [0,tfinal], with dt=tfinal/nsamples, using RK4.
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  orbit_observer observer(&veq, &gc);
  std::size_t nsamples; command_line("nsamples", 512) >> nsamples;
  boost::numeric::odeint::runge_kutta4<guiding_centre::state>
      integration_algorithm;
  boost::numeric::odeint::integrate_const(
      integration_algorithm, odeint_adapter(&gc),
      initial_state, 0.0, tfinal, tfinal/nsamples, observer);

  return 0;
}
