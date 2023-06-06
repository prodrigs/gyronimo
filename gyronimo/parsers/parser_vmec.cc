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

// @parser_vmec.cc, this file is part of ::gyronimo::

#include <gyronimo/core/error.hh>
#include <gyronimo/core/linspace.hh>
#include <gyronimo/parsers/parser_vmec.hh>

namespace gyronimo {

//! Reads and parses a `VMEC` netcdf ouput file by name.
parser_vmec::parser_vmec(const std::string& filename) {
  try {
    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
    get_data(dataFile, "Aminor_p", Aminor_p_);
    get_data(dataFile, "Rmajor_p", Rmajor_p_);
    get_data(dataFile, "aspect", aspect_);
    get_data(dataFile, "b0", b0_);
    get_data(dataFile, "bdotgradv", bdotgradv_);
    get_data(dataFile, "beta_vol", beta_vol_);
    get_data(dataFile, "betapol", beta_pol_);
    get_data(dataFile, "betator", beta_tor_);
    get_data(dataFile, "betatotal", beta_total_);
    get_data(dataFile, "betaxis", beta_axis_);
    get_data(dataFile, "buco", buco_);
    get_data(dataFile, "bvco", bvco_);
    get_data(dataFile, "chi", chi_);
    get_data(dataFile, "iotaf", iotaf_);
    get_data(dataFile, "iotas", iotas_);
    get_data(dataFile, "jcuru", jcuru_);
    get_data(dataFile, "jcurv", jcurv_);
    get_data(dataFile, "jdotb", jdotb_);
    get_data(dataFile, "lasym__logical__", is_axisymmetric_);
    get_data(dataFile, "mass", mass_);
    get_data(dataFile, "mnmax", mnmax_);
    get_data(dataFile, "mnmax_nyq", mnmax_nyq_);
    get_data(dataFile, "mpol", mpol_);
    get_data(dataFile, "nfp", nfp_);
    get_data(dataFile, "ns", ns_);
    get_data(dataFile, "ntor", ntor_);
    get_data(dataFile, "phi", phi_);
    get_data(dataFile, "phipf", phipf_);
    get_data(dataFile, "phips", phips_);
    get_data(dataFile, "pres", pres_);
    get_data(dataFile, "presf", presf_);
    get_data(dataFile, "q_factor", q_factor_);
    get_data(dataFile, "raxis_cc", raxis_cc_);
    get_data(dataFile, "rbtor", rbtor_);
    get_data(dataFile, "rbtor0", rbtor0_);
    get_data(dataFile, "rmax_surf", rmax_surf_);
    get_data(dataFile, "rmin_surf", rmin_surf_);
    get_data(dataFile, "signgs", signgs_);
    get_data(dataFile, "version_", version_);
    get_data(dataFile, "volume_p", volume_p_);
    get_data(dataFile, "xm", xm_);
    get_data(dataFile, "xm_nyq", xm_nyq_);
    get_data(dataFile, "xn", xn_);
    get_data(dataFile, "xn_nyq", xn_nyq_);
    get_data(dataFile, "zaxis_cs", zaxis_cs_);
    get_data(dataFile, "zmax_surf", zmax_surf_);
    get_data_2d(dataFile, "bmnc", bmnc_);
    get_data_2d(dataFile, "bsubsmns", bsubsmns_);
    get_data_2d(dataFile, "bsubumnc", bsubumnc_);
    get_data_2d(dataFile, "bsubvmnc", bsubvmnc_);
    get_data_2d(dataFile, "bsupumnc", bsupumnc_);
    get_data_2d(dataFile, "bsupvmnc", bsupvmnc_);
    get_data_2d(dataFile, "gmnc", gmnc_);
    get_data_2d(dataFile, "lmns", lmns_);
    get_data_2d(dataFile, "rmnc", rmnc_);
    get_data_2d(dataFile, "zmns", zmns_);
  } catch (netCDF::exceptions::NcException& e) {
    std::cout << e.what() << std::endl;
  }
  radius_ = linspace<narray_type>(0.0, 1.0, ns_);
  double ds_half_cell = 0.5 / (ns_ - 1);
  radius_half_cell_ =
      linspace<narray_type>(ds_half_cell, 1.0 - ds_half_cell, ns_ - 1);
}
void parser_vmec::get_data(
    const netCDF::NcFile& nc, const std::string& var, bool& out) {
  auto data = nc.getVar(var);
  if (data.isNull() || data.getDimCount() != 0)
    error(__func__, __FILE__, __LINE__, "Bad data in " + var, 1);
  int dataIn;
  data.getVar(&dataIn);
  out = (dataIn % 2 ? true : false);
}
void parser_vmec::get_data(
    const netCDF::NcFile& nc, const std::string& var, int& out) {
  auto data = nc.getVar(var);
  if (data.isNull() || data.getDimCount() != 0)
    error(__func__, __FILE__, __LINE__, "Bad data in " + var, 1);
  data.getVar(&out);
}
void parser_vmec::get_data(
    const netCDF::NcFile& nc, const std::string& var, size_t& out) {
  auto data = nc.getVar(var);
  if (data.isNull() || data.getDimCount() != 0)
    error(__func__, __FILE__, __LINE__, "Bad data in " + var, 1);
  int dataIn;
  data.getVar(&dataIn);
  out = (size_t)dataIn;
}
void parser_vmec::get_data(
    const netCDF::NcFile& nc, const std::string& var, double& out) {
  auto data = nc.getVar(var);
  if (data.isNull() || data.getDimCount() != 0)
    error(__func__, __FILE__, __LINE__, "Bad data in " + var, 1);
  data.getVar(&out);
}
void parser_vmec::get_data(
    const netCDF::NcFile& nc, const std::string& var, narray_type& out) {
  auto data = nc.getVar(var);
  if (data.isNull() || data.getDimCount() != 1)
    error(__func__, __FILE__, __LINE__, "Bad data in " + var, 1);
  auto dims = data.getDims();
  size_t size = dims[0].getSize();
  out.resize(size);
  data.getVar(&out[0]);
}
void parser_vmec::get_data_2d(
    const netCDF::NcFile& nc, const std::string& var, narray_type& out) {
  auto data = nc.getVar(var);
  if (data.isNull() || data.getDimCount() != 2)
    error(__func__, __FILE__, __LINE__, "Bad data in " + var, 1);
  auto dims = data.getDims();
  size_t size0 = dims[0].getSize();
  size_t size1 = dims[1].getSize();
  out.resize(size0 * size1);
  data.getVar(&out[0]);
}

}  // end namespace gyronimo.
