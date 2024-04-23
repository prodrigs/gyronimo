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

// @vmecdump.cc, this file is part of ::gyronimo::

// Command-line tool to extract info from `VMEC` output files.
// External dependencies:
// - [argh](https://github.com/adishavit/argh), a minimalist argument handler.
// - [GSL](https://www.gnu.org/software/gsl), the GNU Scientific Library.
// - [boost](https://www.boost.org), the boost library.
// - [netcdf-c++4] (https://github.com/Unidata/netcdf-cxx4.git).

#include <gyronimo/core/linspace.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>
#include <gyronimo/metrics/metric_vmec.hh>
#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/version.hh>

#include <argh.h>
#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

using namespace gyronimo;

void print_help() {
  std::cout << "vmecdump, powered by ::gyronimo::v" << version_major << "."
            << version_minor << "." << version_patch
            << " (git-commit:" << git_commit_hash << ").\n";
  std::string help_message =
      "usage: vmecdump [options] vmec_netcdf_file\n"
      "reads a vmec output file, prints required information to stdout.\n"
      "options:\n"
      "  -info  Prints general info about the equilibrium.\n"
      "  -prof  Prints the radial grid, iota, and pressure profiles.\n"
      "  -rphiz Reads u v w triplets from stdin, prints R phi Z to stdout.\n"
      "  -surface [options] [scalar-field, scalar-field,...]\n"
      "         Reads a u sequence from stdin, prints required scalar fields\n"
      "         along the corresponding flux surface, ordered as below:\n"
      "         -r, -z, -phi\n"
      "                Any of the r/z/phi coordinates.\n"
      "         -jac   Metric jacobian.\n"
      "         -b     field magnitude, by IR3field::magnitude().\n"
      "         Available options:\n"
      "         -python\n"
      "                Output in python array format (default table).\n"
      "         -nzeta=n, -ntheta=n\n"
      "                Sets surface sampling rate (default 75).\n";
  std::cout << help_message;
  std::exit(0);
}

void print_info(const parser_vmec& vmec) {
  std::cout << "axisymmetric: " << (vmec.is_axisymmetric() ? "yes\n" : "no\n");
  std::cout << "fieldperiods: " << vmec.nfp() << '\n';
  std::cout << "     nradial: " << vmec.ns() << '\n';
  std::cout << "        ntor: " << vmec.ntor() << '\n';
  std::cout << "        mpol: " << vmec.mpol() << '\n';
  std::cout << "      signgs: " << vmec.signgs() << '\n';
  std::cout << "       B_mag: " << vmec.B0() << " [T]" << '\n';
  std::cout << "       R_mag: " << vmec.R0() << " [m]" << '\n';
  std::cout << "       F_mag: " << vmec.rbtor0() << " [m.T]" << '\n';
  std::cout << "       R_geo: " << vmec.rmajor() << " [m]" << '\n';
  std::cout << "     a_minor: " << vmec.aminor() << " [m]" << '\n';
  std::cout << "      volume: " << vmec.volume() << " [m^3]" << '\n';
}

void print_profiles(const parser_vmec& vmec) {
  for (size_t i = 0; i < vmec.ns(); i++)
    std::cout << vmec.sgrid()[i] << " " << vmec.iotaf()[i] << " "
              << vmec.pres()[i] << '\n';
}

void print_rphiz(
    const parser_vmec& vmap, const interpolator1d_factory* ifactory) {
  dblock_adapter u_range(vmap.sgrid());
  interpolator1d** Rmnc = new interpolator1d*[vmap.xm().size()];
  interpolator1d** Zmns = new interpolator1d*[vmap.xm().size()];
  for (size_t i = 0; i < vmap.xm().size(); i++) {
    std::slice u_slice(i, u_range.size(), vmap.xm().size());
    std::valarray<double> rmnc_i = (vmap.rmnc())[u_slice];
    Rmnc[i] = ifactory->interpolate_data(u_range, dblock_adapter(rmnc_i));
    std::valarray<double> zmns_i = (vmap.zmns())[u_slice];
    Zmns[i] = ifactory->interpolate_data(u_range, dblock_adapter(zmns_i));
  };
  double u, v, w;
  while (std::cin >> u >> v >> w) {
    double R = 0.0, Z = 0.0;
    for (size_t i = 0; auto m : vmap.xm()) {
      double n = vmap.xn()[i];
      R += (*Rmnc[i])(u)*std::cos(m * w - n * v);
      Z += (*Zmns[i])(u)*std::sin(m * w - n * v);
      i++;
    }
    std::cout << R << " " << v << " " << Z << "\n";
  };
}

void print_surface(
    const parser_vmec& vmap, const interpolator1d_factory* ifactory,
    const argh::parser& command_line) {
  morphism_vmec morph(&vmap, ifactory);
  metric_vmec g(&morph);
  equilibrium_vmec veq(&g, ifactory);
  size_t nzeta, ntheta;
  command_line("nzeta", 75) >> nzeta;
  command_line("ntheta", 75) >> ntheta;
  auto zeta_range =
      linspace<parser_vmec::narray_type>(0.0, 2 * std::numbers::pi, nzeta);
  auto theta_range =
      linspace<parser_vmec::narray_type>(0.0, 2 * std::numbers::pi, ntheta);
  double u;
  while (std::cin >> u) {
    if (u <= 0.0 || u > 1.0) continue;  // ignores invalid s values.
    for (double w : theta_range) {
      for (double v : zeta_range) {
        const IR3 q = {u, v, w};
        auto [R, z] = morph.get_rz(q);
        double phi = v, x = R * std::cos(phi), y = R * std::sin(phi);
        if (!command_line["python"])
          std::cout << x << " " << y << " " << z << " ";
        if (command_line["r"]) std::cout << R << " ";
        if (command_line["z"]) std::cout << z << " ";
        if (command_line["phi"]) std::cout << phi << " ";
        if (command_line["jac"]) std::cout << g.jacobian(q) << " ";
        if (command_line["b"]) std::cout << veq.magnitude(q, 0.0) << " ";
        if (!command_line["python"]) std::cout << '\n';
        else if (v < 2 * std::numbers::pi) std::cout << " ";
      }
      std::cout << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {
  auto command_line = argh::parser(argv);
  if (command_line[{"h", "help"}]) print_help();
  if (!command_line(1)) {  // the 1st non-option argument is the vmec file.
    std::cout << "vmecdump: no vmec_netcdf file provided; -h for help.\n";
    std::exit(1);
  }
  parser_vmec vmec(command_line[1]);
  cubic_gsl_factory ifactory;
  if (command_line["info"]) print_info(vmec);
  if (command_line["prof"]) print_profiles(vmec);
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  if (command_line["rphiz"]) print_rphiz(vmec, &ifactory);
  if (command_line["surface"]) print_surface(vmec, &ifactory, command_line);
  return 0;
}
