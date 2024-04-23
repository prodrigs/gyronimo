// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues.

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

// @heltrace.cc, this file is part of ::gyronimo::

// Command-line tool to print guiding-centre orbits in `HELENA` equilibria.
// External dependencies:
// - [argh](https://github.com/adishavit/argh), a minimalist argument handler.
// - [GSL](https://www.gnu.org/software/gsl), the GNU Scientific Library.
// - [boost](https://www.boost.org), the boost library.

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/linspace.hh>
#include <gyronimo/dynamics/guiding_centre.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>
#include <gyronimo/fields/equilibrium_helena.hh>
#include <gyronimo/interpolators/bicubic_gsl.hh>
#include <gyronimo/parsers/parser_helena.hh>
#include <gyronimo/version.hh>

#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <argh.h>
#include <cmath>
#include <iostream>

void print_help() {
  std::cout << "heltrace, powered by ::gyronimo::v" << gyronimo::version_major
            << "." << gyronimo::version_minor << "." << gyronimo::version_patch
            << " (git-commit:" << gyronimo::git_commit_hash << ").\n";
  std::string help_message =
      "usage: heltrace [options] helena_output_file\n"
      "reads an helena output file, prints guiding-centre orbit to stdout.\n"
      "options:\n"
      "  -rhom= Axis density (in m_proton*1e19, default 1).\n"
      "  -mass= Particle mass (in m_proton, default 1).\n"
      "  -pphi= Pphi value (in eV.s, default 1).\n"
      "  -charge=\n"
      "         Particle charge (in q_proton, default 1).\n"
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
  orbit_observer(
      double zstar, double vstar, const gyronimo::IR3field_c1* e,
      const gyronimo::guiding_centre* g)
      : zstar_(zstar), vstar_(vstar), eq_pointer_(e), gc_pointer_(g) {};
  void operator()(const gyronimo::guiding_centre::state& s, double t) {
    gyronimo::IR3 x = gc_pointer_->get_position(s);
    double v_parallel = gc_pointer_->get_vpp(s);
    double bphi = eq_pointer_->covariant_versor(x, t)[gyronimo::IR3::w];
    double flux = x[gyronimo::IR3::u] * x[gyronimo::IR3::u];
    std::cout << t << " " << x[gyronimo::IR3::u] << " " << x[gyronimo::IR3::v]
              << " " << x[gyronimo::IR3::w] << " " << v_parallel << " "
              << -zstar_ * flux + vstar_ * v_parallel * bphi << " "
              << gc_pointer_->energy_perpendicular(s, t) << " "
              << gc_pointer_->energy_parallel(s) << '\n';
  };
 private:
  double zstar_, vstar_;
  const gyronimo::IR3field_c1* eq_pointer_;
  const gyronimo::guiding_centre* gc_pointer_;
};

// Finds the radial position s = \sqrt{\Psi/\Psi_b} at the midplane.
/*
   The canonical toroidal momentum is
     P_phi = q_s A_\phi + m_s v_\parallel B_\phi/B,
   where
     \mathbf{B} = \nabla \phi \times \nabla \Psi + B_\phi \nabla \phi,
     \mathbf{B} = \nabla \times \mathbf{A} => A_\phi = - \Psi,
   and all variables are in SI units, except \Psi in Wb/rad. Writting P_\phi in
   eV.s and the energy E in eV, with e the electron charge, one gets:
     P_phi = - (q_s/e) \Psi +
         \sigma_\parallel b_\phi \sqrt{2 m_s (E/e) (1 - \Lambda B/B_0)}
   Function arguments:
   pphi: \P_\phi in eV.s;
   zstar: (q_s/e) \Psi_b, with \Psi_b the boundary flux, q_s the species charge;
   vdagger: \sqrt{2 m_s (E/e)}, with E in eV, m_s the species mass;
   lambda: \Lambda;
   vpp_sign: the name says it all;
   heq: reference to an HELENA equilibrium object;
*/
double get_initial_radial_position(
    double pphi, double zstar, double vdagger, double lambda, double vpp_sign,
    const gyronimo::equilibrium_helena& heq) {
  auto pphi_functional = [zstar, vdagger, lambda, vpp_sign, &heq](double s) {
    gyronimo::IR3 pos = {s, 0.0, 0.0};  // Assuming {s,chi,phi} coordinates!
    double Btilde = heq.magnitude(pos, 0.0);
    double bphi = heq.covariant_versor(pos, 0.0)[gyronimo::IR3::w];
    return -zstar * s * s +
        vpp_sign * vdagger * bphi * std::sqrt(1.0 - lambda * Btilde);
  };
  auto orbit = [&pphi_functional, pphi](double s) {
    return pphi_functional(s) - pphi;
  };
  auto s_grid = gyronimo::linspace<std::valarray<double>>(0.0, 1.0, 1024);
  auto bracketing_iterator = std::adjacent_find(
      std::begin(s_grid), std::end(s_grid),
      [&orbit](double x, double y) { return orbit(x) * orbit(y) < 0.0; });
  if (bracketing_iterator == std::end(s_grid)) {
    std::cout << "# orbit not crossing the low-field side midplane.\n"
              << "# try pphi in the range [" << pphi_functional(0.0) << ":"
              << pphi_functional(1.0) << "]\n";
    std::exit(1);
  }
  auto root_interval = boost::math::tools::bisect(
      [&orbit](double s) { return orbit(s); }, *bracketing_iterator,
      *(bracketing_iterator + 1),
      [](double a, double b) { return std::abs(b - a) < 1.0e-09; });
  return 0.5 * (root_interval.first + root_interval.second);
}

int main(int argc, char* argv[]) {
  auto command_line = argh::parser(argv);
  if (command_line[{"h", "help"}]) print_help();
  if (!command_line(1)) {  // the 1st non-option argument is the mapping file.
    std::cout << "heltrace: no helena mapping file provided; -h for help.\n";
    std::exit(1);
  }
  gyronimo::parser_helena hmap(command_line[1]);
  gyronimo::bicubic_gsl_factory ifactory(
      false, (hmap.is_symmetric() ? 0 : 9), (hmap.is_symmetric() ? 9 : 0));
  gyronimo::morphism_helena morph(&hmap, &ifactory);
  gyronimo::metric_helena g(&morph, &ifactory);
  gyronimo::equilibrium_helena heq(&g, &ifactory);

  double pphi, mass, rhom, charge, energy, lambda, tfinal;
  command_line("pphi", 1.0) >> pphi;  // pphi in eV.s.
  command_line("mass", 1.0) >> mass;  // m_proton units.
  command_line("rhom", 1.0) >> rhom;  // density in m_proton*1e19.
  command_line("charge", 1.0) >> charge;  // q_electron units.
  command_line("energy", 1.0) >> energy;  // energy in eV.
  command_line("lambda", 1.0) >> lambda;  // lambda is signed!
  command_line("tfinal", 1.0) >> tfinal;
  double vpp_sign = std::copysign(1.0, lambda);  // lambda carries vpp sign.
  lambda = std::abs(lambda);  // once vpp sign is stored, lambda turns unsigned.

  // Computes normalisation constants:
  double Valfven = heq.B0() /
      std::sqrt(gyronimo::codata::mu0 *
                (rhom * gyronimo::codata::m_proton * 1.e+19));
  double Ualfven = 0.5 * gyronimo::codata::m_proton * mass * Valfven * Valfven;
  double energySI = energy * gyronimo::codata::e;
  double Lref = heq.R0();

  // Prints output header:
  std::cout << "heltrace, powered by ::gyronimo::v" << gyronimo::version_major
            << "." << gyronimo::version_minor << "." << gyronimo::version_patch
            << " (git-commit:" << gyronimo::git_commit_hash << ").\n";
  std::cout << "# args: ";
  for (int i = 1; i < argc; i++) std::cout << argv[i] << " ";
  std::cout << std::endl;
  std::cout << "# l_ref = " << Lref << " [m];";
  std::cout << " v_alfven = " << Valfven << " [m/s];";
  std::cout << " u_alfven = " << Ualfven << " [J];";
  std::cout << " energy = " << energySI << " [J]." << "\n";
  std::cout << "# vars: t s chi phi vpar Pphi/e Eperp/Ealfven Epar/Ealfven\n";

  // Builds the guiding_centre object:
  gyronimo::guiding_centre gc(
      Lref, Valfven, charge / mass, lambda * energySI / Ualfven, &heq, nullptr);

  // Computes the initial conditions from the supplied constants of motion:
  double zstar = charge * g.parser()->cpsurf() * heq.B0() * heq.R0() * heq.R0();
  double vstar =
      Valfven * mass * gyronimo::codata::m_proton / gyronimo::codata::e;
  double vdagger = vstar * std::sqrt(energySI / Ualfven);
  double initial_radial_position =
      get_initial_radial_position(pphi, zstar, vdagger, lambda, vpp_sign, heq);
  gyronimo::guiding_centre::state initial_state = gc.generate_state(
      {initial_radial_position, 0.0, 0.0}, energySI / Ualfven,
      (vpp_sign > 0 ? gyronimo::guiding_centre::plus :
                      gyronimo::guiding_centre::minus),
      0);

  // integrates for t in [0,tfinal], with dt=tfinal/nsamples, using RK4.
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  orbit_observer observer(zstar, vstar, &heq, &gc);
  std::size_t nsamples;
  command_line("samples", 512) >> nsamples;
  boost::numeric::odeint::runge_kutta4<gyronimo::guiding_centre::state>
      integration_algorithm;
  boost::numeric::odeint::integrate_const(
      integration_algorithm, gyronimo::odeint_adapter(&gc), initial_state, 0.0,
      tfinal, tfinal / nsamples, observer);

  return 0;
}
