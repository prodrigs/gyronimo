// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Jorge Ferreira and Paulo Rodrigues.

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
// - [boost](https://www.gnu.org/software/gsl), the boost library.
// - [netcdf-c++4] (https://github.com/Unidata/netcdf-cxx4.git).

#include <cmath>
#include <argh.h>
#include <numbers>
#include <iostream>
#include <gyronimo/version.hh>
#include <gyronimo/core/linspace.hh>
#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/metrics/metric_vmec.hh>
#include <gyronimo/fields/equilibrium_vmec.hh>
#include <gyronimo/interpolators/cubic_gsl.hh>

void print_help() {
  std::cout << "vmecdump, powered by ::gyronimo:: v"
      << gyronimo::version_major << "." << gyronimo::version_minor << ".\n";
  std::cout <<
      "usage: vmecdump [options] vmec_netcdf_file\n";
  std::cout <<
      "reads a VMEC netcdf file and prints the required orbit to stdout.\n";
  std::cout << "options:\n";
  std::cout << "  -rphiz          Reads a {s, u, v} sequence from stdin and prints {R, phi, Z} to stdout.\n";
  std::cout << "  -surface        Reads a {s} sequence from stdin and prints surface{s}(theta,zeta) to stdout for the scalar field:\n";
  std::cout << "      -Bnorm             - magnitude of magnetic field (default option)\n";
  std::cout << "      -R                 - R coordinate\n";
  std::cout << "      -phi               - phi coordinate\n";
  std::cout << "      -Z                 - Z coordinate\n";
  std::cout << "      -jacobian          - jacobian of the metric\n";
  std::cout << "      -jacobian_vmec     - jacobian of the metric (from vmec)\n";
  std::cout << "  -radius         Prints radial grid.\n";
  std::cout << "  -iota           Prints 'iota' profile.\n";
  std::cout << "  -q              Prints 'q' safety factor profile.\n";
  std::cout << "  -pressure       Prints 'pressure' profile.\n";
  std::cout << "  -python_format  Outputs profiles in python array format.\n";
  std::exit(0);
}

enum Scalar { R_cyl, phi_cyl, Z_cyl, Bnorm, Bnorm_vmec, jacobian, jacobian_vmec };
typedef std::valarray<double> narray_type;
enum class Out_format {table, python};

int main(int argc, char* argv[]) {
  void debug(const gyronimo::parser_vmec&);
  void print_rphiz(const gyronimo::parser_vmec&, const gyronimo::interpolator1d_factory*);
  void print_surface(const Scalar, const gyronimo::parser_vmec&, 
                     const gyronimo::interpolator1d_factory*,
                     const int ntheta, const int nzeta);
  void print_surface(const Scalar, const gyronimo::parser_vmec&, const gyronimo::interpolator1d_factory*);
  void printNarray (const gyronimo::parser_vmec::narray_type&, 
                    const int, 
                    const std::string&,
                    const Out_format&);

  auto command_line = argh::parser(argv);
  if (command_line[{"h", "help"}]) print_help();
  if (!command_line(1)) {  // the 1st non-option argument is the VMEC file.
    std::cout << "vmecdump: no VMEC netcdf file provided; -h for help.\n";
    std::exit(1);
  }
  gyronimo::parser_vmec vmec(command_line[1]);
  gyronimo::cubic_gsl_factory ifactory;
 
 // Reads parameters from the command line:
  if (command_line["rphiz"]) {
    print_rphiz(vmec, &ifactory);
    std::exit(0);
  }
  if (command_line["surface"]) {
    int ntheta; command_line("ntheta", 128) >> ntheta;
    int nzeta; command_line("nzeta", 128) >> nzeta;
    Scalar scalar = Bnorm;
    if (command_line["Bnorm_vmec"]) scalar = Bnorm_vmec;
    if (command_line["R"]) scalar = R_cyl;
    if (command_line["Z"]) scalar = Z_cyl;
    if (command_line["jacobian"]) scalar = jacobian;
    if (command_line["jacobian_vmec"]) scalar = jacobian_vmec;

    print_surface(scalar, vmec, &ifactory, ntheta, nzeta);
    std::exit(0);
  }
  double radius; command_line("radius", 1.0) >> radius;  // vmec normalized radius

  Out_format print_format = ((command_line["python_format"]) ? Out_format::python : Out_format::table);

  std::cout << "\n# VMEC equilibrium (" << command_line[1] << ")" << ::std::endl;
  std::cout << "# is axisymmetric?: " << (vmec.is_axisymmetric() ? "No" : "Yes") << std::endl;
  std::cout << "# radial grid size: " << vmec.ns() << std::endl;
  std::cout << "#       grid: rmin = " << vmec.rmin_surf() << ", rmax = " << vmec.rmax_surf() << std::endl;
  std::cout << "#         ns: " << vmec.ns() << ", ntor: " << vmec.ntor() << ", mpol: " << vmec.mpol() << std::endl;
  std::cout << "#        nfp: " << vmec.nfp() << std::endl;
  std::cout << "#     signgs: " << vmec.signgs() << std::endl;
  std::cout << "#      B_mag: " << vmec.B_0() << " [T]" << std::endl;
  std::cout << "#      R_mag: " << vmec.R_0() << " [m]" << std::endl;
  std::cout << "#      F_mag: " << vmec.rbtor0() << " [m.T]" << std::endl;
  std::cout << "#    a_minor: " << vmec.Aminor() << " [m]" << std::endl;
  std::cout << "#       Rgeo: " << vmec.Rmajor() << " [m]" << std::endl;
  std::cout << "# tot.volume: " << vmec.volume() << " [m^3]" << std::endl;
  std::cout << "#_________________________________________________________\n" << std::endl;
   
  if (command_line["radius"]) printNarray(vmec.radius(), 4, "radius", print_format);
  if (command_line["iota"]) printNarray(vmec.iotaf(), 4, "iota", print_format);
  if (command_line["q"]) printNarray(vmec.q(), 4, "q", print_format);
  if (command_line["pressure"]) printNarray(vmec.pres(), 4, "pressure", print_format);
  if (command_line[{"debug"}]) debug(vmec);

  return 0;
}

void print_rphiz(const gyronimo::parser_vmec& vmap, const gyronimo::interpolator1d_factory *ifactory) {
  gyronimo::dblock_adapter s_range(vmap.radius());
  // gyronimo::dblock_adapter u_range(linspace<narray_type>(0.0,std::numbers::pi , ns_));
  gyronimo::interpolator1d **Rmnc = new gyronimo::interpolator1d* [vmap.xm().size()];
  gyronimo::interpolator1d **Zmns = new gyronimo::interpolator1d* [vmap.xm().size()];
  for(size_t i=0; i<vmap.xm().size(); i++) {
    // std::cout << i <<"\n";
      std::slice s_cut (i, s_range.size(), vmap.xm().size());
      std::valarray<double> rmnc_i = (vmap.rmnc())[s_cut];
      Rmnc[i] = ifactory->interpolate_data(s_range, gyronimo::dblock_adapter(rmnc_i));
      std::valarray<double> zmns_i = (vmap.zmns())[s_cut];
      Zmns[i] = ifactory->interpolate_data(s_range, gyronimo::dblock_adapter(zmns_i));
      // for (auto ss: s_range) std::cout << (*Rmnc_[i])(ss) << ", ";
      // std::cout << "----------------------" <<std::endl;
  };
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  double u, v, w;
  while (std::cin >> u >> v >> w ) {
    // std::cout << s << " " << u << " "<< v << "\n";
    double R = 0.0, Z = 0.0;
    for (size_t i = 0; auto m : vmap.xm()) {  
      double n = vmap.xn()[i];
      R+= (*Rmnc[i])(u) * std::cos( m*w - n*v ); 
      Z+= (*Zmns[i])(u) * std::sin( m*w - n*v );
      // std::cout << i <<" " << R << " " << v << " " << Z << "\n";
      i++;
    }
    std::cout << R << " " << v << " " << Z << "\n";
  };
}

// print surface in theta, zeta
void print_surface(const Scalar scalar,
                   const gyronimo::parser_vmec& vmap, 
                   const gyronimo::interpolator1d_factory *ifactory,
                   const int ntheta, const int nzeta) {

  gyronimo::metric_vmec g(&vmap, ifactory);
  gyronimo::equilibrium_vmec veq(&g, ifactory);
  gyronimo::dblock_adapter s_range(vmap.radius());
  gyronimo::interpolator1d **Rmnc = new gyronimo::interpolator1d* [vmap.xm().size()];
  gyronimo::interpolator1d **Zmns = new gyronimo::interpolator1d* [vmap.xm().size()];
  for(size_t i=0; i<vmap.xm().size(); i++) {
      std::slice s_cut (i, s_range.size(), vmap.xm().size());
      std::valarray<double> rmnc_i = (vmap.rmnc())[s_cut];
      Rmnc[i] = ifactory->interpolate_data(s_range, gyronimo::dblock_adapter(rmnc_i));
      std::valarray<double> zmns_i = (vmap.zmns())[s_cut];
      Zmns[i] = ifactory->interpolate_data(s_range, gyronimo::dblock_adapter(zmns_i));
  };
  double u;
  gyronimo::dblock_adapter theta_range(gyronimo::linspace<narray_type>(0.0, 2*std::numbers::pi, ntheta));
  gyronimo::dblock_adapter zeta_range(gyronimo::linspace<narray_type>(0.0, 2*std::numbers::pi, nzeta));
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  while (std::cin >> u ) {
    for(double w: theta_range) {
      for(auto j=0; double v: zeta_range) {
        double R = 0.0, Z = 0.0;
        for (size_t i = 0; i < vmap.xm().size(); i++) {  
          double m = vmap.xm()[i]; double n = vmap.xn()[i];
          R+= (*Rmnc[i])(u) * std::cos( m*w - n*v ); 
          Z+= (*Zmns[i])(u) * std::sin( m*w - n*v );
        }
        const gyronimo::IR3 p = {u,v,w};  // {s, zeta, theta}
        switch(scalar)
        {
          case Bnorm: std::cout << veq.magnitude(p, 0.0); break;
          case Bnorm_vmec: std::cout << veq.magnitude_vmec(p, 0.0); break;
          case R_cyl: std::cout << R ; break;
          case Z_cyl: std::cout << Z ; break;
          case jacobian: std::cout << g.jacobian(p); break;
          case jacobian_vmec: std::cout << g.jacobian_vmec(p); break;
        }
        if(++j < nzeta) std::cout << " ";
      }
      std::cout << std::endl;
    };
  };
}
void print_surface(const Scalar scalar, const gyronimo::parser_vmec& vmap, const gyronimo::interpolator1d_factory *ifactory) {
  print_surface(scalar, vmap, ifactory, 128, 128);
}
//print printNarray with n columns 
void printNarray (const gyronimo::parser_vmec::narray_type& va, 
                  const int ncols, 
                  const std::string&  label,
                  const Out_format& format)
{
    std::cout.precision(8);
    std::cout.setf(std::ios::scientific);
    if (format == Out_format::python) {
    // python format
      std::cout << "\n" << label << " = [ " << std::endl;
      for (int i=0; i<(va.size()/ncols); i++) {
          for (int j=0; j<ncols; j++) {
              std::cout << va[i*ncols+j];
              if ((i*ncols+j)<(va.size()-1)) std::cout << ", ";
          }
          std::cout << std::endl;
      }  
      for (int i=(va.size()-va.size()%ncols); i<va.size(); i++) {
        std::cout << va[i];
        if (i<(va.size()-1)) std::cout << ", ";
      }
      std::cout <<  "]\n" << std::endl;
    }
    // table format
    else {
      std::cout << "\n# " << label << " : " << std::endl;
      for (int i=0; i<(va.size()/ncols); i++) {
          for (int j=0; j<ncols; j++) {
              std::cout << va[i*ncols+j] << " ";
          }  
          std::cout << std::endl;
      }
      for (int i=(va.size()-va.size()%ncols); i<va.size(); i++)
      {
        std::cout << va[i] << " ";
      }
      std::cout << std::endl;
    };
}

void debug(const gyronimo::parser_vmec& vmec ) {
  auto print_format = Out_format::python;
  double u, v, w;
  std::cin >> u >> v >> w;
  gyronimo::IR3 x {u, v, w};
  gyronimo::cubic_gsl_factory ifactory;
  gyronimo::metric_vmec g(&vmec, &ifactory);
  auto metric_at_x = g(x);
  auto dmetric_at_x = g.del(x);
  auto jacobian_at_x = g.jacobian(x);
  auto J_at_x = g.jacobian_vmec(x);
  gyronimo::equilibrium_vmec veq(&g, &ifactory);
  auto B_at_x = veq.contravariant(x, 0.0);
  auto norm_b = veq.magnitude(x, 0.0);
  auto d_i_B_at_x = veq.del_contravariant(x, 0.0);
  std::cout << "x = ("
        << x[gyronimo::IR3::u] << ", "
        << x[gyronimo::IR3::v] << ", "
        << x[gyronimo::IR3::w] << ") "
        << std::endl;
  std::cout << "g = ["
        << metric_at_x[gyronimo::SM3::uu] << ", "
        << metric_at_x[gyronimo::SM3::uv] << ", "
        << metric_at_x[gyronimo::SM3::uw] << ", "
        << metric_at_x[gyronimo::SM3::vv] << ", "
        << metric_at_x[gyronimo::SM3::vw] << ", "
        << metric_at_x[gyronimo::SM3::ww] << "] "
        << std::endl;
  std::cout << "dg = ["
      << dmetric_at_x[gyronimo::dSM3::uuu] << ", "
      << dmetric_at_x[gyronimo::dSM3::uuv] << ", "
      << dmetric_at_x[gyronimo::dSM3::uuw] << ", "
      << dmetric_at_x[gyronimo::dSM3::uvu] << ", "
      << dmetric_at_x[gyronimo::dSM3::uvv] << ", "
      << dmetric_at_x[gyronimo::dSM3::uvw] << ", "
      << dmetric_at_x[gyronimo::dSM3::uwu] << ", "
      << dmetric_at_x[gyronimo::dSM3::uwv] << ", "
      << dmetric_at_x[gyronimo::dSM3::uww] << ", "
      << dmetric_at_x[gyronimo::dSM3::vvu] << ", "
      << dmetric_at_x[gyronimo::dSM3::vvv] << ", "
      << dmetric_at_x[gyronimo::dSM3::vvw] << ", "
      << dmetric_at_x[gyronimo::dSM3::vwu] << ", "
      << dmetric_at_x[gyronimo::dSM3::vwv] << ", "
      << dmetric_at_x[gyronimo::dSM3::vww] << ", "
      << dmetric_at_x[gyronimo::dSM3::wwu] << ", "
      << dmetric_at_x[gyronimo::dSM3::wwv] << ", "
      << dmetric_at_x[gyronimo::dSM3::www] << "] "
      << std::endl;
  std::cout << "J = " << jacobian_at_x << std::endl;
  std::cout << "J (vmec)= " << J_at_x << std::endl;
  std::cout << "B = ("
        << B_at_x[gyronimo::IR3::u] << ", "
        << B_at_x[gyronimo::IR3::v] << ", "
        << B_at_x[gyronimo::IR3::w] << ") "
        << std::endl;
  std::cout << "norm_B = " << norm_b << std::endl;
  std::cout << "dB = ("
        << B_at_x[gyronimo::dIR3::uu] << ", "
        << B_at_x[gyronimo::dIR3::uv] << ", "
        << B_at_x[gyronimo::dIR3::uw] << ", "
        << B_at_x[gyronimo::dIR3::vu] << ", "
        << B_at_x[gyronimo::dIR3::vv] << ", "
        << B_at_x[gyronimo::dIR3::vw] << ", "
        << B_at_x[gyronimo::dIR3::wu] << ", "
        << B_at_x[gyronimo::dIR3::wv] << ", "        
        << B_at_x[gyronimo::dIR3::ww] << ") "
        << std::endl;
  std::cout << "cpsurf = " << vmec.cpsurf() << std::endl;
  printNarray(vmec.xm(), 4, "xm", print_format);
  printNarray(vmec.xn(), 4, "xn", print_format);
  printNarray(vmec.xm_nyq(), 4, "xm_nyq", print_format);
  printNarray(vmec.xn_nyq(), 4, "xn_nyq", print_format);
  printNarray(vmec.rmnc(), 4, "rmnc", print_format);
  printNarray(vmec.zmns(), 4, "zmns", print_format);
  printNarray(vmec.rmnc(), 4, "rmnc", print_format);
  printNarray(vmec.gmnc(), 4, "gmnc", print_format);
  printNarray(vmec.bmnc(), 4, "bmnc", print_format);
  printNarray(vmec.lmns(), 4, "lmns", print_format);
}
