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

// @parser_vmec.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_PARSER_VMEC
#define GYRONIMO_PARSER_VMEC

#include <netcdf>
#include <string>
#include <valarray>

namespace gyronimo {

//! Parsing object for `VMEC` output files.
/*!
    Reads and parses the output file produced by `VMEC`, a Grad-Shafanov
    equilibrium code [S.P. Hirshman et al., Phys. Fluids **26**, 3353 (1983)]
    distributed with the stellarator optimisation package
    [STELLOPT](https://princetonuniversity.github.io/STELLOPT/VMEC.html). The
    output is a netcdf file (usually named `wout`) and the full description of
    stored fields can be found by running `ncdump -h wout_filename`. Parsed data
    is accessed via public member functions, some of them returning
    `parser_vmec::narray_type` arrays.
*/
class parser_vmec {
 public:
  typedef std::valarray<double> narray_type;
  parser_vmec(const std::string& filename);
  ~parser_vmec() {};

  int signgs() const { return signgs_; };
  bool is_axisymmetric() const { return is_axisymmetric_; };
  size_t mnmax() const { return mnmax_; };
  size_t mnmax_nyq() const { return mnmax_nyq_; };
  size_t mpol() const { return mpol_; };
  size_t nfp() const { return nfp_; };
  size_t ns() const { return ns_; };
  size_t ntor() const { return ntor_; };
  size_t version() const { return version_; };
  double aminor() const { return Aminor_p_; };
  double B0() const { return b0_; };
  double R0() const { return rbtor0_ / b0_; };
  double rmajor() const { return Rmajor_p_; };
  double aspect() const { return aspect_; };
  double beta_axis() const { return beta_axis_; };
  double beta_pol() const { return beta_pol_; };
  double beta_tor() const { return beta_tor_; };
  double beta_total() const { return beta_total_; };
  double cpsurf() const;
  double rbtor() const { return rbtor_; };
  double rbtor0() const { return rbtor0_; };
  double rmax_surf() const { return rmax_surf_; };
  double rmin_surf() const { return rmin_surf_; };
  double volume() const { return volume_p_; };
  double zmax_surf() const { return zmax_surf_; };
  const narray_type& bdotb() const { return bdotb_; };
  const narray_type& bdotgradv() const { return bdotgradv_; };
  const narray_type& beta_vol() const { return beta_vol_; };
  const narray_type& bmnc() const { return bmnc_; };
  const narray_type& bsubsmns() const { return bsubsmns_; };
  const narray_type& bsubumnc() const { return bsubumnc_; };
  const narray_type& bsubvmnc() const { return bsubvmnc_; };
  const narray_type& bsupumnc() const { return bsupumnc_; };
  const narray_type& bsupvmnc() const { return bsupvmnc_; };
  const narray_type& buco() const { return buco_; };
  const narray_type& bvco() const { return bvco_; };
  const narray_type& chi() const { return chi_; };
  const narray_type& chipf() const { return chipf_; };
  const narray_type& currumnc() const { return currumnc_; };
  const narray_type& currvmnc() const { return currvmnc_; };
  const narray_type& gmnc() const { return gmnc_; };
  const narray_type& iotaf() const { return iotaf_; };
  const narray_type& iotas() const { return iotas_; };
  const narray_type& jcuru() const { return jcuru_; };
  const narray_type& jcurv() const { return jcurv_; };
  const narray_type& jdotb() const { return jdotb_; };
  const narray_type& lmns() const { return lmns_; };
  const narray_type& mass() const { return mass_; };
  const narray_type& phi() const { return phi_; };
  const narray_type& phip() const { return phip_; };
  const narray_type& phipf() const { return phipf_; };
  const narray_type& phips() const { return phips_; };
  const narray_type& pres() const { return pres_; };
  const narray_type& presf() const { return presf_; };
  const narray_type& q() const { return q_factor_; };
  const narray_type& sgrid() const { return sgrid_; };
  const narray_type& sgrid_half_cell() const { return sgrid_half_cell_; };
  const narray_type& raxis_cc() const { return raxis_cc_; };
  const narray_type& rmnc() const { return rmnc_; };
  const narray_type& xm() const { return xm_; };
  const narray_type& xm_nyq() const { return xm_nyq_; };
  const narray_type& xn() const { return xn_; };
  const narray_type& xn_nyq() const { return xn_nyq_; };
  const narray_type& zaxis_cs() const { return zaxis_cs_; };
  const narray_type& zmns() const { return zmns_; };
 private:
  int signgs_;
  bool is_axisymmetric_;
  size_t mnmax_;
  size_t mnmax_nyq_;
  size_t mpol_;
  size_t nfp_;
  size_t ns_;
  size_t ntor_;
  size_t version_;
  double Aminor_p_;
  double Rmajor_p_;
  double aspect_;
  double b0_;
  double beta_axis_;
  double beta_pol_;
  double beta_tor_;
  double beta_total_;
  double rbtor0_;
  double rbtor_;
  double rmax_surf_;
  double rmin_surf_;
  double volume_p_;
  double zmax_surf_;
  narray_type bdotb_;
  narray_type bdotgradv_;
  narray_type beta_;
  narray_type beta_vol_;
  narray_type bmnc_;
  narray_type bsubsmns_;
  narray_type bsubumnc_;
  narray_type bsubvmnc_;
  narray_type bsupumnc_;
  narray_type bsupvmnc_;
  narray_type buco_;
  narray_type bvco_;
  narray_type chi_;
  narray_type chipf_;
  narray_type currumnc_;
  narray_type currvmnc_;
  narray_type gmnc_;
  narray_type iota_;
  narray_type iotaf_;
  narray_type iotas_;
  narray_type jcuru_;
  narray_type jcurv_;
  narray_type jdotb_;
  narray_type lmns_;
  narray_type mass_;
  narray_type phi_;
  narray_type phip_;
  narray_type phipf_;
  narray_type phips_;
  narray_type pres_;
  narray_type presf_;
  narray_type q_factor_;
  narray_type sgrid_;
  narray_type sgrid_half_cell_;
  narray_type raxis_cc_;
  narray_type rmnc_;
  narray_type xm_;
  narray_type xm_nyq_;
  narray_type xn_;
  narray_type xn_nyq_;
  narray_type zaxis_cs_;
  narray_type zmns_;

  void get_data(const netCDF::NcFile&, const std::string&, int&);
  void get_data(const netCDF::NcFile&, const std::string&, bool&);
  void get_data(const netCDF::NcFile&, const std::string&, size_t&);
  void get_data(const netCDF::NcFile&, const std::string&, double&);
  void get_data(const netCDF::NcFile&, const std::string&, narray_type&);
  void get_data_2d(const netCDF::NcFile&, const std::string&, narray_type&);
};

inline double parser_vmec::cpsurf() const {
  return chi_[chi_.size() - 1] / (rbtor0_ * rbtor0_ / b0_);
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_PARSER_VMEC
