// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues, Jorge Ferreira.

// @parser_vmec.hh

#pragma once

#include <string>
#include <iostream>
#include <valarray>
#include <netcdf>
#include <stdio.h>
#include <assert.h> 

namespace gyronimo {
//! Parsing object for `VMEC` mapping files.
/*!
    Reads and parses a NetCDF file produced by `VMEC`, ... 
*/
class parser_vmec {
 public:
  typedef std::valarray<double> narray_type;
  parser_vmec(const std::string& filename);
  ~parser_vmec() {};
  bool is_axisymmetric() const {return is_axisymmetric_;};
  size_t version()       const {return version_;};
  size_t nfp()           const {return nfp_;};
  size_t ns()            const {return ns_;};
  size_t nradius()       const {return nradius_;};
  size_t ntor()          const {return ntor_;};
  size_t mpol()          const {return mpol_;};
  size_t mnmax()         const {return mnmax_;};
  size_t mnmax_nyq()     const {return mnmax_nyq_;};
  int signgs()           const {return signgs_;};
  // VMEC scalars
  double aspect_ratio() const { return aspect_ratio_;};
  double Aminor()       const { return Aminor_p_;};
  double Rmajor()       const { return Rmajor_p_;};
  double volume()       const { return volume_p_;};
  double rmin_surf()    const { return rmin_surf_;};
  double rmax_surf()    const { return rmax_surf_;};
  double zmax_surf()    const { return zmax_surf_;};
  double beta_total()   const { return beta_total_;};
  double beta_pol()     const { return beta_pol_;};
  double beta_tor()     const { return beta_tor_;};
  double beta_axis()    const { return beta_axis_;};
  double B_0()          const { return b0_;};
  double R_0()          const { return rbtor0_/b0_;};
  double cpsurf()       const { return chi_[chi_.size()-1]/(rbtor0_*rbtor0_/b0_);};
  double rbtor0()       const { return rbtor0_;};
  double rbtor()        const { return rbtor_;};
// VMEC radial grids
  const narray_type& radius()        const { return radius_;};
  const narray_type& radius_half()   const { return radius_half_;};
// VMEC profiles
  const narray_type& iotaf()    const { return iotaf_;};
  const narray_type& q()        const { return q_factor_;};
  const narray_type& presf()    const { return presf_;};
  const narray_type& phi()      const { return phi_; };
  const narray_type& phip()     const { return phip_; };
  const narray_type& phipf()    const { return phipf_;};
  const narray_type& chi()      const { return chi_;};
  const narray_type& chipf()    const { return chipf_;};
  const narray_type& mass()     const { return mass_;};
  const narray_type& pres()     const { return pres_;};
  const narray_type& iotas()    const { return iotas_;};
  const narray_type& beta_vol() const { return beta_vol_;};
  const narray_type& phips()    const { return phips_;};
  const narray_type& bvco()     const { return bvco_;};
  const narray_type& buco()     const { return buco_;};
  const narray_type& jcuru()    const { return jcuru_;};
  const narray_type& jcurv()    const { return jcurv_;};
  const narray_type& jdotb()    const { return jdotb_;};
  const narray_type& bdotb()    const { return bdotb_;};
  const narray_type& bdotgradv()const { return bdotgradv_;};
  // VMEC spectral representation
  const narray_type& raxis_cc() const { return raxis_cc_;};
  const narray_type& zaxis_cs() const { return zaxis_cs_;};
  const narray_type& xm()       const { return xm_;};
  const narray_type& xn()       const { return xn_;};
  const narray_type& xm_nyq()   const { return xm_nyq_;};
  const narray_type& xn_nyq()   const { return xn_nyq_;};
  const narray_type& rmnc()     const { return rmnc_;};
  const narray_type& zmns()     const { return zmns_; };
  const narray_type& lmns()     const { return lmns_; };
  const narray_type& gmnc()     const { return gmnc_; };
  const narray_type& bmnc()     const { return bmnc_; };
  const narray_type& bsubumnc() const { return bsubumnc_; };
  const narray_type& bsubvmnc() const { return bsubvmnc_; };
  const narray_type& bsubsmns() const { return bsubsmns_; };
  const narray_type& currumnc() const { return currumnc_; };
  const narray_type& currvmnc() const { return currvmnc_; };
  const narray_type& bsupumnc() const { return bsupumnc_; };
  const narray_type& bsupvmnc() const { return bsupvmnc_; };

  // in VMEC coordinates
  double ns_b_;
  double nfp_b_;
  double aspect_ratio_b_;
  double rmin_b_;
  double rmax_b_;
  double beta_axis_b;

 private:
  bool is_axisymmetric_;
  size_t version_;
  size_t nfp_, ns_;
  size_t nradius_;
  size_t ntor_;
  size_t mpol_;
  size_t mnmax_;
  size_t mnmax_nyq_;
  int signgs_;
  double aspect_ratio_;
  double Aminor_p_;
  double Rmajor_p_;
  double volume_p_;
  double rmin_surf_;
  double rmax_surf_;
  double zmax_surf_;
  double beta_total_;
  double beta_pol_;
  double beta_tor_;
  double beta_axis_;
  double b0_;
  double rbtor0_;
  double rbtor_;
// VMEC profiles
  narray_type radius_;
  narray_type radius_half_;
  narray_type iota_;
  narray_type pres_;
  narray_type beta_;
  narray_type phip_; 
  narray_type phi_;
  narray_type bvco_, buco_;
  narray_type iotaf_;
  narray_type q_factor_;
  narray_type presf_;
  narray_type phipf_;
  narray_type chi_;
  narray_type chipf_;
  narray_type mass_;
  narray_type iotas_;
  narray_type beta_vol_;
  narray_type phips_;
  narray_type jcuru_;
  narray_type jcurv_;
  narray_type jdotb_;
  narray_type bdotb_;
  narray_type bdotgradv_;
  narray_type raxis_cc_;
  narray_type zaxis_cs_;
  narray_type xm_;
  narray_type xn_;
  narray_type xm_nyq_;
  narray_type xn_nyq_;
  narray_type rmnc_;
  narray_type zmns_;
  narray_type lmns_;
  narray_type gmnc_;
  narray_type bmnc_;
  narray_type bsubumnc_;
  narray_type bsubvmnc_;
  narray_type bsubsmns_;
  narray_type currumnc_;
  narray_type currvmnc_;
  narray_type bsupumnc_;
  narray_type bsupvmnc_;

  void getData(const netCDF::NcFile&, const std::string&, bool&);
  void getData(const netCDF::NcFile&, const std::string&, int&);
  void getData(const netCDF::NcFile&, const std::string&, size_t&);
  void getData(const netCDF::NcFile&, const std::string&, double&);
  void getData(const netCDF::NcFile&, const std::string&, narray_type&);
  void getData2D(const netCDF::NcFile&, const std::string&,  narray_type&);
}; 

} // end namespace gyronimo.
