// :: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Jorge Ferreira and Paulo Rodrigues.

// :: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// :: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::.  If not, see <https://www.gnu.org/licenses/>.

// @vmecdump.cc, this file is part of ::

// Command-line tool to extract info from `VMEC` output files.
// External dependencies:
// - [argh](https://github.com/adishavit/argh), a minimalist argument handler.
// - [GSL](https://www.gnu.org/software/gsl), the GNU Scientific Library.
// - [boost](https://www.boost.org), the boost library.
// - [netcdf-c++4] (https://github.com/Unidata/netcdf-cxx4.git).

#include <cmath>
#include <argh.h>
#include <string>
#include <numbers>
#include <iostream>
#include <gyronimo/version.hh>
#include <gyronimo/core/linspace.hh>
#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/metrics/metric_vmec.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>

void print_help() {
  std::cout << "vmecdump, powered by ::v"
      << gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
  std::string help_message =
      "usage: vmecdump [options] vmec_netcdf_file\n"
      "reads a vmec output file, prints the required information to stdout.\n"
      "options:\n"
      "  -info     Prints general info about the equilibrium.\n"
      "  -prof     Prints the radial grid, iota, and pressure profiles.\n"
      "  -rphiz    Reads u v w lines from stdin, prints R phi Z to stdout.\n"
      "  -surface  Reads a u sequence from stdin, prints the required scalar\n"
      "            field along the corresponding flux surface.\n"
      "    -b      B magnitude, by IR3field::magnitude() (default).\n"
      "    -b_vmec    B magnitude, via vmec output.\n"
      "    -r/-z/-phi Any of the R/Z/phi coordinates.\n"
      "    -jac       Metric jacobian, by metric_covariant::jacobian().\n"
      "    -jac_vmec  Native metric jacobian, via vmec output.\n"
      "    -python    Output in python array format (default table).\n"
      "    -nzeta=/-ntheta= Sets surface sampling rate (default 75).\n";
  std::cout << help_message;
  std::exit(0);
}

using namespace gyronimo;
enum format_t {table, python};
enum scalar_t {r_cyl, phi_cyl, z_cyl, b_gymo, b_vmec, jac_gymo, jac_vmec};

int main(int argc, char* argv[]) {
  void print_info(const parser_vmec&);
  void print_profiles(const parser_vmec&);
  void print_rphiz(const parser_vmec&, const interpolator1d_factory*);
  void print_surface(
      const parser_vmec&, const interpolator1d_factory*,
      const int, const int, const scalar_t, const format_t);
  auto command_line = argh::parser(argv);
  if(command_line[{"h", "help"}]) print_help();
  if(!command_line(1)) {  // the 1st non-option argument is the vmec file.
    std::cout << "vmecdump: no vmec_netcdf file provided; -h for help.\n";
    std::exit(1);
  }
  parser_vmec vmec(command_line[1]);
  cubic_gsl_factory ifactory;
  if(command_line["info"]) print_info(vmec);
  if(command_line["prof"]) print_profiles(vmec);
  if(command_line["rphiz"]) print_rphiz(vmec, &ifactory);
  if(command_line["surface"]) {
    int nzeta; command_line("nzeta", 75) >> nzeta;
    int ntheta; command_line("ntheta", 75) >> ntheta;
    format_t format = (command_line["python"] ? python : table);
    scalar_t scalar = b_gymo;
    if(command_line["r"]) scalar = r_cyl;
    if(command_line["z"]) scalar = z_cyl;
    if(command_line["phi"]) scalar = phi_cyl;
    if(command_line["jac"]) scalar = jac_gymo;
    if(command_line["b_vmec"]) scalar = b_vmec;
    if(command_line["jac_vmec"]) scalar = jac_vmec;
    print_surface(vmec, &ifactory, ntheta, nzeta, scalar, format);
  }
  return 0;
}
void print_info(const parser_vmec& vmec) {
  std::cout << "axisymmetric: " << (vmec.is_axisymmetric() ? "yes\n" : "no\n");
  std::cout << "fieldperiods: " << vmec.nfp() << "\n";
  std::cout << "     nradial: " << vmec.ns() << "\n";
  std::cout << "        ntor: " << vmec.ntor() << "\n";
  std::cout << "        mpol: " << vmec.mpol() << "\n";
  std::cout << "      signgs: " << vmec.signgs() << "\n";
  std::cout << "       B_mag: " << vmec.B_0() << " [T]" << "\n";
  std::cout << "       R_mag: " << vmec.R_0() << " [m]" << "\n";
  std::cout << "       F_mag: " << vmec.rbtor0() << " [m.T]" << "\n";
  std::cout << "       R_geo: " << vmec.Rmajor() << " [m]" << "\n";
  std::cout << "     a_minor: " << vmec.Aminor() << " [m]" << "\n";
  std::cout << "      volume: " << vmec.volume() << " [m^3]" << "\n";
}
void print_profiles(const parser_vmec& vmec) {
  for(size_t i = 0;i < vmec.ns();i++)
      std::cout << vmec.radius()[i] << " "
          << vmec.iotaf()[i] << " " << vmec.pres()[i] << "\n";
}
void print_rphiz(
    const parser_vmec& vmap, const interpolator1d_factory *ifactory) {
  dblock_adapter u_range(vmap.radius());
  interpolator1d **Rmnc = new interpolator1d* [vmap.xm().size()];
  interpolator1d **Zmns = new interpolator1d* [vmap.xm().size()];
  for(size_t i=0; i<vmap.xm().size(); i++) {
      std::slice u_slice (i, u_range.size(), vmap.xm().size());
      std::valarray<double> rmnc_i = (vmap.rmnc())[u_slice];
      Rmnc[i] = ifactory->interpolate_data(u_range, dblock_adapter(rmnc_i));
      std::valarray<double> zmns_i = (vmap.zmns())[u_slice];
      Zmns[i] = ifactory->interpolate_data(u_range, dblock_adapter(zmns_i));
  };
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  double u, v, w;
  while(std::cin >> u >> v >> w ) {
    double R = 0.0, Z = 0.0;
    for(size_t i = 0; auto m : vmap.xm()) {  
      double n = vmap.xn()[i];
      R += (*Rmnc[i])(u)*std::cos(m*w - n*v); 
      Z += (*Zmns[i])(u)*std::sin(m*w - n*v);
      i++;
    }
    std::cout << R << " " << v << " " << Z << "\n";
  };
}
void print_surface(
    const parser_vmec& vmap,
    const interpolator1d_factory *ifactory,
    const int ntheta, const int nzeta,
    const scalar_t scalar, const format_t format) {
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  metric_vmec g(&vmap, ifactory);
  equilibrium_vmec veq(&g, ifactory);
  auto zeta_range =
      linspace<parser_vmec::narray_type>(0.0, 2*std::numbers::pi, nzeta);
  auto theta_range =
      linspace<parser_vmec::narray_type>(0.0, 2*std::numbers::pi, ntheta);
  double u; while(std::cin >> u ) {
    for(double w : theta_range) {
      for(double v : zeta_range) {
        const IR3 p = {u, v, w};
        IR3 c = veq.metric()->transform2cylindrical(p);
        double R = c[IR3::u], phi = c[IR3::v], z = c[IR3::w];
        double x = R*std::cos(phi), y = R*std::sin(phi);
        switch(format) {
          case table: std::cout << x << " " << y << " " << z << " ";break;
          case python: break;
        }
        switch(scalar) {
          case b_gymo: std::cout << veq.magnitude(p, 0.0); break;
          case b_vmec: std::cout << veq.magnitude_vmec(p, 0.0); break;
          case r_cyl: std::cout << R ; break;
          case z_cyl: std::cout << z ; break;
          case phi_cyl: std::cout << phi ; break;
          case jac_gymo: std::cout << g.jacobian(p); break;
          case jac_vmec: std::cout << g.jacobian_vmec(p); break;
        }
        switch(format) {
          case table: std::cout << std::endl; break;
          case python: if(v < 2*std::numbers::pi) std::cout << " "; break;
        }
      }
      std::cout << std::endl;
    }
  }
}
