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

// @heldump.cc, this file is part of ::gyronimo::

// Command line tool to extract info from `HELENA` output files.
// External dependencies:
// - [argh](https://github.com/adishavit/argh), a minimalist argument handler;
// - [GSL](https://www.gnu.org/software/gsl), the GNU Scientific Library;

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/dblock.hh>
#include <gyronimo/core/linspace.hh>
#include <gyronimo/core/transpose.hh>
#include <gyronimo/interpolators/bicubic_gsl.hh>
#include <gyronimo/parsers/parser_helena.hh>
#include <gyronimo/version.hh>

#include <argh.h>
#include <cmath>
#include <iostream>
#include <numbers>

void print_help() {
  std::cout << "heldump, powered by ::gyronimo::v" << gyronimo::version_major
            << "." << gyronimo::version_minor << "." << gyronimo::version_patch
            << " (git-commit:" << gyronimo::git_commit_hash << ").\n";
  std::string help_message =
      "usage: heldump [options] helena_output_file\n"
      "reads an helena output file, prints required information to stdout.\n"
      "options:\n"
      "  -info  Prints general info.\n"
      "  -prof  Prints 1d profiles: s, pressure, curj, q, dq/ds.\n"
      "  -rz    Reads a {s, chi} sequence from stdin, prints RZ to stdout.\n"
      "  -levels\n"
      "         Reads a s sequence from stdin and prints the corresponding\n"
      "         surfaces (-nchi=N poloidal samples, default 128) to stdout.\n";
  std::cout << help_message;
  std::exit(0);
}

void print_info(const gyronimo::parser_helena& hmap) {
  std::cout << "symmetric: " << (hmap.is_symmetric() ? "yes" : "no") << "\n";
  std::cout << "npoloidal: " << hmap.nchi() << "\n";
  std::cout << "  nradial: " << hmap.npsi() << "\n";
  std::cout << "  eps_geo: " << hmap.eps() << " aka epsilon=a_geo/R_geo.\n";
  std::cout << "  eps_mag: " << hmap.radius() << " aka radius=a_geo/R_mag.\n";
  std::cout << "   cpsurf: " << hmap.cpsurf() << " psi_B/(B_mag*R_mag^2).\n";
  std::cout << "    P_mag: "
            << hmap.p0()[0] * hmap.bmag() * hmap.bmag() / gyronimo::codata::mu0
            << " [Pa].\n";
  std::cout << "    B_mag: " << hmap.bmag() << " [T].\n";
  std::cout << "    R_mag: " << hmap.rmag() << " [m].\n";
  std::cout << "    R_geo: " << hmap.rgeo() << " [m].\n";
}

void print_prof(const gyronimo::parser_helena& hmap) {
  for (size_t i = 0; i < hmap.npsi(); i++)
    std::cout << hmap.s()[i] << " " << hmap.p0()[i] / hmap.p0()[0] << " "
              << hmap.curj()[i] << " " << hmap.q()[i] << " " << hmap.dqs()[i]
              << "\n";
}

void print_rz(const gyronimo::parser_helena& hmap) {
  gyronimo::dblock_adapter s_range(hmap.s());
  gyronimo::dblock_adapter chi_range(hmap.chi());
  gyronimo::bicubic_gsl x(
      s_range, chi_range, gyronimo::dblock_adapter(hmap.x()), false,
      (hmap.is_symmetric() ? 0 : 9), (hmap.is_symmetric() ? 9 : 0));
  gyronimo::bicubic_gsl y(
      s_range, chi_range, gyronimo::dblock_adapter(hmap.y()), false,
      (hmap.is_symmetric() ? 0 : 9), (hmap.is_symmetric() ? 9 : 0));
  double s, chi;
  while (std::cin >> s >> chi) {
    chi -= 2 * std::numbers::pi * std::floor(chi / (2 * std::numbers::pi));
    double y_value;
    if (hmap.is_symmetric() && chi > std::numbers::pi) {
      chi = 2 * std::numbers::pi - chi;
      y_value = -y(s, chi);
    } else y_value = y(s, chi);
    std::cout << hmap.rgeo() * (1.0 + hmap.eps() * x(s, chi)) << " "
              << hmap.rgeo() * hmap.eps() * y_value << "\n";
  }
}

void print_levels(
    const gyronimo::parser_helena& hmap, const argh::parser& command_line) {
  gyronimo::dblock_adapter s_range(hmap.s()), chi_range(hmap.chi());
  gyronimo::bicubic_gsl x(
      s_range, chi_range, gyronimo::dblock_adapter(hmap.x()), false,
      (hmap.is_symmetric() ? 0 : 9), (hmap.is_symmetric() ? 9 : 0));
  gyronimo::bicubic_gsl y(
      s_range, chi_range, gyronimo::dblock_adapter(hmap.y()), false,
      (hmap.is_symmetric() ? 0 : 9), (hmap.is_symmetric() ? 9 : 0));
  size_t nchi;
  command_line("nchi", 128) >> nchi;
  double delta_chi = 2 * std::numbers::pi / nchi;
  auto chi_array = gyronimo::linspace<gyronimo::parser_helena::narray_type>(
      0.0, nchi * delta_chi, nchi);
  double s;
  while (std::cin >> s) {
    if (s <= 0.0 || s > 1.0) continue;  // ignores invalid s values.
    for (double chi : chi_array) {
      chi -= 2 * std::numbers::pi * std::floor(chi / (2 * std::numbers::pi));
      double y_value;
      if (hmap.is_symmetric() && chi > std::numbers::pi) {
        chi = 2 * std::numbers::pi - chi;
        y_value = -y(s, chi);
      } else y_value = y(s, chi);
      std::cout << hmap.rgeo() * (1.0 + hmap.eps() * x(s, chi)) << " "
                << hmap.rgeo() * hmap.eps() * y_value << "\n";
    }
    std::cout << "\n";
  }
}

int main(int argc, char* argv[]) {
  auto command_line = argh::parser(argv);
  if (command_line[{"h", "help"}]) print_help();
  if (!command_line(1)) {  // the 1st non-option argument is the mapping file.
    std::cout << "heldump: no helena mapping file provided; -h for help.\n";
    std::exit(1);
  }
  gyronimo::parser_helena hmap(command_line[1]);
  if (command_line["rz"]) print_rz(hmap);
  if (command_line["info"]) print_info(hmap);
  if (command_line["prof"]) print_prof(hmap);
  if (command_line["levels"]) print_levels(hmap, command_line);
  return 0;
}
