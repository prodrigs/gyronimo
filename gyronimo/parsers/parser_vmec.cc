// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues, Jorge Ferreira.

// @parser_vmec.cc

#include <algorithm>
#include <vector>
#include <numbers>
#include <fstream>
#include <random>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/io_streaming.hh>
#include <gyronimo/parsers/parser_vmec.hh>
#include <gyronimo/core/linspace.hh>

namespace gyronimo {
  typedef std::valarray<double> narray_type;
  //! Reads and parses a VMEC netcdf ouput file by name.
  parser_vmec::parser_vmec(const std::string& filename) {
    try {
      netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
      getData(dataFile, "version_", version_);
      getData(dataFile, "lasym__logical__", is_axisymmetric_);
      getData(dataFile, "ns", ns_);
      nradius_ = ns_; 
      radius_ = linspace<narray_type>(0.0, 1.0, ns_);
      double ds_half_cell = 0.5/(ns_-1);
      radius_half_ = linspace<narray_type>(ds_half_cell, 1.0-ds_half_cell, ns_-1);
      // VMEC scalars
      getData(dataFile, "nfp", nfp_);
      getData(dataFile, "ntor", ntor_);
      getData(dataFile, "mpol", mpol_);
      getData(dataFile, "signgs", signgs_);
      getData(dataFile, "mnmax", mnmax_);
      getData(dataFile, "mnmax_nyq", mpol_);
      getData(dataFile, "aspect", mnmax_nyq_);
      getData(dataFile, "Aminor_p", Aminor_p_);
      getData(dataFile, "Rmajor_p", Rmajor_p_);
      getData(dataFile, "volume_p", volume_p_);
      getData(dataFile, "rmin_surf", rmin_surf_);
      getData(dataFile, "rmax_surf", rmax_surf_);
      getData(dataFile, "zmax_surf", zmax_surf_);
      getData(dataFile, "betatotal", beta_total_);
      getData(dataFile, "betapol", beta_pol_);
      getData(dataFile, "betator", beta_tor_);
      getData(dataFile, "betaxis", beta_axis_);
      getData(dataFile, "b0", b0_);
      getData(dataFile, "rbtor0", rbtor0_);
      getData(dataFile, "rbtor", rbtor_);
      // VMEC profiles
      getData(dataFile, "iotaf", iotaf_);
      getData(dataFile, "q_factor", q_factor_);
      getData(dataFile, "presf", presf_);
      getData(dataFile, "phi", phi_);
      getData(dataFile, "phipf", phipf_);
      getData(dataFile, "chi", chi_); //  not present in vacuum solutions
      getData(dataFile, "mass", mass_);
      getData(dataFile, "pres", pres_);
      getData(dataFile, "iotas", iotas_);
      getData(dataFile, "beta_vol", beta_vol_);
      getData(dataFile, "phips", phips_);
      getData(dataFile, "bvco", bvco_);
      getData(dataFile, "buco", buco_);
      getData(dataFile, "jcuru", jcuru_);
      getData(dataFile, "jcurv", jcurv_);
      getData(dataFile, "jdotb", jdotb_);
      getData(dataFile, "bdotgradv", bdotgradv_);
      // profiles in ntor
      getData(dataFile, "raxis_cc", raxis_cc_);
      getData(dataFile, "zaxis_cs", zaxis_cs_);
      // spectral related data
      getData(dataFile, "xm", xm_);
      getData(dataFile, "xn", xn_);
      getData(dataFile, "xm_nyq", xm_nyq_);
      getData(dataFile, "xn_nyq", xn_nyq_);
      getData2D(dataFile, "rmnc", rmnc_);
      getData2D(dataFile, "zmns", zmns_);
      getData2D(dataFile, "lmns", lmns_);
      getData2D(dataFile, "gmnc", gmnc_);
      getData2D(dataFile, "bmnc", bmnc_);
      getData2D(dataFile, "bsubumnc", bsubumnc_);
      getData2D(dataFile, "bsubvmnc", bsubvmnc_);
      getData2D(dataFile, "bsubsmns", bsubsmns_);
      getData2D(dataFile, "bsupumnc", bsupumnc_);
      getData2D(dataFile, "bsupvmnc", bsupvmnc_);
    } catch (netCDF::exceptions::NcException &e) {
      std::cout << e.what() << std::endl;
    }
  }
  void parser_vmec::getData(const netCDF::NcFile& nc, const std::string& var, bool& out){
      auto data = nc.getVar(var); 
      if(data.isNull()) error(__func__, __FILE__, __LINE__,
                        ("Null data for var: "+var).c_str() , 1);
      assert( data.getDimCount() == 0 );
      int dataIn; data.getVar(&dataIn);
      out = ( dataIn % 2 ? true : false );
  }
  void parser_vmec::getData(const netCDF::NcFile& nc, const std::string& var, int& out){
      auto data = nc.getVar(var); 
      if(data.isNull()) error(__func__, __FILE__, __LINE__,
                  ("Null data for var: "+var).c_str() , 1);assert( data.getDimCount() == 0 );
      data.getVar(&out);
  }
  void parser_vmec::getData(const netCDF::NcFile& nc, const std::string& var, size_t& out){
      auto data = nc.getVar(var); 
      if(data.isNull()) error(__func__, __FILE__, __LINE__,
                  ("Null data for var: "+var).c_str() , 1);
      assert( data.getDimCount() == 0 );
      int dataIn; data.getVar(&dataIn);
      out = (size_t) dataIn;
  }
  void parser_vmec::getData(const netCDF::NcFile& nc, const std::string& var, double& out){
      auto data = nc.getVar(var); 
      if(data.isNull()) error(__func__, __FILE__, __LINE__,
                  ("Null data for var: "+var).c_str() , 1);
      assert( data.getDimCount() == 0 );
      data.getVar(&out);
  }
  void parser_vmec::getData(const netCDF::NcFile& nc, const std::string& var,  narray_type& out){
      auto data = nc.getVar(var); 
      if(data.isNull()) error(__func__, __FILE__, __LINE__,
                  ("Null data for var: "+var).c_str() , 1);
      assert( data.getDimCount() == 1 );
      auto dims = data.getDims(); size_t size  = dims[0].getSize();
      out.resize(size);
      data.getVar(&out[0]);
  }
  void parser_vmec::getData2D(const netCDF::NcFile& nc, const std::string& var,  narray_type& out){
      auto data = nc.getVar(var); 
      if(data.isNull()) error(__func__, __FILE__, __LINE__,
                  ("Null data for var: "+var).c_str() , 1);
      assert( data.getDimCount() == 2 );
      auto dims = data.getDims(); 
      size_t size0  = dims[0].getSize(); size_t size1  = dims[1].getSize();
      out.resize(size0*size1);
      data.getVar(&out[0]);
  }
} // end namespace gyronimo.
