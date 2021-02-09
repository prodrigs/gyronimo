// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @heldump.cc
// Command line tool to extract info from `HELENA` output files.
// External dependencies:
// - [argh](https://github.com/adishavit/argh), a minimalist argument handler;
// - [GSL](https://www.gnu.org/software/gsl), the GNU Scientific Library;

#include <cmath>
#include <numbers>
#include <iostream>
#include <argh/argh.h>
#include <gyronimo/core/codata.hh>
#include <gyronimo/core/dblock.hh>
#include <gyronimo/core/version.hh>
#include <gyronimo/core/linspace.hh>
#include <gyronimo/core/transpose.hh>
#include <gyronimo/parsers/parser_helena.hh>
#include <gyronimo/interpolators/bicubic_gsl.hh>

void print_help() {
  std::cout << "heldump, powered by ::gyronimo:: v"
      << gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
  std::cout << "usage: heldump [options] hmap\n";
  std::cout <<
  "reads an hmap file produced by HELENA and prints required info to stdout.\n";
  std::cout << "options:\n";
  std::cout <<
  "    -info   Prints general info.\n" <<
  "    -prof   Prints 1d profiles: s, pressure, curj, q, dq/ds.\n" <<
  "    -rz     Reads a {s, chi} sequence from stdin, prints RZ to stdout.\n" <<
  "    -xy     Reads a {s, chi} sequence from stdin, prints XY to stdout.\n" <<
  "    -levels Reads a s sequence from stdin and prints the corresponding\n" <<
  "            surfaces (-nchi=N poloidal samples, default 128) to stdout.\n";
  std::exit(0);
}

int main(int argc, char *argv[]) {
  void print_info(const gyronimo::parser_helena&);
  void print_prof(const gyronimo::parser_helena&);
  void print_xy_rz(const gyronimo::parser_helena&, const argh::parser&);
  void print_levels(const gyronimo::parser_helena&, const argh::parser&);
  auto command_line = argh::parser(argv);
  if (command_line[{"h", "help"}]) print_help();
  if (!command_line(1)) {  // the 1st non-option argument is the mapping file.
    std::cout << "heldump: no helena mapping file provided; -h for help.\n";
    std::exit(1);
  }
  gyronimo::parser_helena hmap(command_line[1]);
  if (command_line["info"]) print_info(hmap);
  if (command_line["prof"]) print_prof(hmap);
  if (command_line["levels"]) print_levels(hmap, command_line);
  if (command_line["rz"] || command_line["xy"]) print_xy_rz(hmap, command_line);
  return 0;
}
void print_info(const gyronimo::parser_helena& hmap) {
  std::cout << "npoloidal: " << hmap.nchi() << "\n";
  std::cout << "nradial:   " << hmap.npsi() << "\n";
  std::cout << "eps_geo: " << hmap.eps() << " aka epsilon = a_geo/R_geo.\n";
  std::cout << "eps_mag: " << hmap.radius() << " aka radius = a_geo/R_mag.\n";
  std::cout << "cpsurf:  " << hmap.cpsurf() << " psi_B/(B_mag*R_mag^2).\n";
  std::cout << "P_mag:   " <<
      hmap.p0()[0]*hmap.bmag()*hmap.bmag()/gyronimo::codata::mu0
          << " (Pa).\n";
  std::cout << "R_mag:   " << hmap.rmag() << " (m).\n";
  std::cout << "B_mag:   " << hmap.bmag() << " (T).\n";
}
void print_prof(const gyronimo::parser_helena& hmap) {
  for (size_t i = 0; i < hmap.npsi(); i++)
    std::cout << hmap.s()[i] << " " << hmap.p0()[i] << " "
        << hmap.curj()[i] << " " << hmap.q()[i] << " "
            << hmap.dqs()[i] << "\n";
}
void print_xy_rz(
    const gyronimo::parser_helena& hmap, const argh::parser& command_line) {
  gyronimo::dblock_adapter s_range(hmap.s());
  gyronimo::dblock_adapter chi_range(hmap.chi());
  gyronimo::bicubic_gsl x(s_range, chi_range,
      gyronimo::dblock_adapter(gyronimo::transpose(hmap.x(), hmap.nchi())));
  gyronimo::bicubic_gsl y(s_range, chi_range,
      gyronimo::dblock_adapter(gyronimo::transpose(hmap.y(), hmap.nchi())));
  double s, chi;
  while (std::cin >> s >> chi) {
    chi -= 2*std::numbers::pi*std::floor(chi/(2*std::numbers::pi));
    double y_val;
    if (hmap.is_symmetric() && chi > std::numbers::pi) {
      chi = 2*std::numbers::pi - chi;
      y_val = -y(s, chi);
    } else
      y_val = y(s, chi);
    if (command_line["xy"]) std::cout << x(s, chi) << " " << y_val << "\n";
    if (command_line["rz"])
        std::cout << hmap.rmag()*(1.0 + hmap.eps()*x(s, chi)) << " "
          << hmap.rmag()*hmap.eps()*y_val << "\n";
  }
}
void print_levels(
    const gyronimo::parser_helena& hmap, const argh::parser& command_line) {
  gyronimo::dblock_adapter s_range(hmap.s()), chi_range(hmap.chi());
  gyronimo::bicubic_gsl x(s_range, chi_range,
      gyronimo::dblock_adapter(gyronimo::transpose(hmap.x(), hmap.nchi())));
  gyronimo::bicubic_gsl y(s_range, chi_range,
      gyronimo::dblock_adapter(gyronimo::transpose(hmap.y(), hmap.nchi())));
  size_t nchi;
  command_line("nchi", 128) >> nchi;  // defaults to 128.
  double delta_chi = 2*std::numbers::pi/nchi;
  auto chi_array =
      gyronimo::linspace<gyronimo::parser_helena::narray_type>(
          0.0, nchi*delta_chi, nchi);
  double s;
  while (std::cin >> s) {
    if (s <=0.0 || s > 1.0) continue;  // ignores invalid s values.
    for (double chi : chi_array) {
      chi -= 2*std::numbers::pi*std::floor(chi/(2*std::numbers::pi));
      double y_val;
      if (hmap.is_symmetric() && chi > std::numbers::pi) {
        chi = 2*std::numbers::pi - chi;
        y_val = -y(s, chi);
      } else y_val = y(s, chi);
      std::cout << x(s, chi) << " " << y_val << "\n";
    }
    std::cout << "\n";
  }
}
