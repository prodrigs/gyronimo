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

// @equilibrium_a_field.cc, this file is part of ::gyronimo::


//to compile: g++ -c -std=c++20 -O2 equilibrium_a_field.cc -I ../gyronimo 

#include <gyronimo/core/dblock.hh>
#include <gyronimo/fields/equilibrium_vmec_a.hh>

#include <cmath>
#include <numeric>

namespace gyronimo {

equilibrium_a_field::equilibrium_a_field(
    const metric_vmec* g, const interpolator1d_factory* ifactory, bool normalised)
    : IR3field_c1(std::abs(g->my_parser()->B0()), 1.0, g),
      metric_(g), parser_(g->my_parser()), harmonics_(parser_->mnmax_nyq()),
      m_(parser_->xm_nyq()), n_(parser_->xn_nyq()), index_(harmonics_),
      btheta_mn_(parser_->mnmax_nyq()), bzeta_mn_(parser_->mnmax_nyq()), g_mn_(parser_->mnmax_nyq()), normalised_(normalised) {

  std::iota(index_.begin(), index_.end(), 0);

  this->build_interpolator_array(
      bzeta_mn_,  parser_->bsupvmnc() / this->m_factor(), ifactory, harmonics_);
  this->build_interpolator_array(
      btheta_mn_, parser_->bsupumnc() / this->m_factor(), ifactory, harmonics_);
  this->build_interpolator_array(
      g_mn_, parser_->gmnc(), ifactory, harmonics_);


  this->build_harmonics(mc_, nc_, m_, n_, csupumnc, csupvmnc);

  index_.resize(mc_.size());
  std::iota(index_.begin(), index_.end(), 0);


  czeta_mn_.resize(mc_.size());
  this->build_interpolator_array(
      czeta_mn_, csupvmnc, ifactory, mc_.size());
  ctheta_mn_.resize(mc_.size());
  this->build_interpolator_array(
      ctheta_mn_, csupumnc, ifactory, mc_.size());


  this->build_integrals(ctheta_integ, czeta_integ);

  dblock_adapter sgrid(parser_->sgrid());
  ctheta_00_ = std::unique_ptr<interpolator1d>( ifactory->interpolate_data( sgrid, dblock_adapter(ctheta_integ) ) );
  czeta_00_ = std::unique_ptr<interpolator1d>( ifactory->interpolate_data( sgrid, dblock_adapter(czeta_integ) ) );
}



IR3 equilibrium_a_field::covariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];

  const auto& cis_mn = equilibrium_a_field::cached_cis_new(theta, zeta);

  double as = 0;
  double azeta = 0;
  double atheta = 0;

  for(size_t i=0; i<mc_.size(); i++) {
    double cos_mn = std::real(cis_mn[i]);
    double sin_mn = std::imag(cis_mn[i]);
    if ( mc_[i]!=0 && nc_[i]==0 ){
        as = as + (*czeta_mn_[i])(s)/mc_[i]*sin_mn; 
    }
    else if ( mc_[i]==0 && nc_[i]!=0 ){
        as = as + (*ctheta_mn_[i])(s)/nc_[i]*sin_mn;
    }
    else if ( mc_[i]==0 && nc_[i]==0 ){
        azeta = (*ctheta_00_)(s);
        atheta = -(*czeta_00_)(s);
    }
    else {
        as = as + ( (*czeta_mn_[i])(s)/mc_[i]*sin_mn + (*ctheta_mn_[i])(s)/nc_[i]*sin_mn )/2;
    }
  }
  return {as/2, azeta/2, atheta/2};
}



IR3 equilibrium_a_field::contravariant(const IR3& position, double time) const {
  IR3 A = this->covariant(position, time);
  return metric_->to_contravariant(A, position);    
}



dIR3 equilibrium_a_field::del_covariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];

  const auto& cis_mn = equilibrium_a_field::cached_cis_new(theta, zeta);
  
  double das_ds = 0;
  double das_dzeta = 0;
  double das_dtheta = 0;
  double dazeta_ds = 0;
  double datheta_ds = 0;

  for(size_t i=0; i<mc_.size(); i++) {
    double cos_mn = std::real(cis_mn[i]);
    double sin_mn = std::imag(cis_mn[i]);
    if ( mc_[i]!=0 && nc_[i]==0 ){
        das_ds = das_ds + (*czeta_mn_[i]).derivative(s)/mc_[i]*sin_mn;
        das_dzeta = das_dzeta - (*czeta_mn_[i])(s)/mc_[i]*nc_[i]*cos_mn;
        das_dtheta = das_dtheta + (*czeta_mn_[i])(s)*cos_mn;
    }
    else if ( mc_[i]==0 && nc_[i]!=0 ){
        das_ds = das_ds + (*ctheta_mn_[i]).derivative(s)/nc_[i]*sin_mn;
        das_dzeta = das_dzeta - (*ctheta_mn_[i])(s)*cos_mn;
        das_dtheta = das_dtheta + (*ctheta_mn_[i])(s)/nc_[i]*mc_[i]*cos_mn;
    }
    else if ( mc_[i]==0 && nc_[i]==0 ){
        dazeta_ds = (*ctheta_mn_[i])(s);
        datheta_ds = -(*czeta_mn_[i])(s);
    }
    else {
        das_ds = das_ds + ( (*czeta_mn_[i]).derivative(s)/mc_[i]*sin_mn + (*ctheta_mn_[i]).derivative(s)/nc_[i]*sin_mn )/2;
        das_dzeta = das_dzeta - ( (*czeta_mn_[i])(s)/mc_[i]*nc_[i]*cos_mn + (*ctheta_mn_[i])(s)*cos_mn )/2;
        das_dtheta = das_dtheta + ( (*czeta_mn_[i])(s)*cos_mn + (*ctheta_mn_[i])(s)/nc_[i]*mc_[i]*cos_mn )/2;
    }
  }
  return {das_ds/2, das_dzeta/2, das_dtheta/2,
          dazeta_ds/2, 0, 0,
          datheta_ds/2, 0, 0};
}



dIR3 equilibrium_a_field::del_contravariant(const IR3& position, double time) const {
    const metric_covariant *g = this->metric();

    dIR3 c1 = contraction<second>(
      (*g).del_inverse(position), this->covariant(position, time));
    dIR3 c2 = contraction<first>(
      this->del_covariant(position, time), (*g).inverse(position));
    return {
      c1[dIR3::uu] + c2[dIR3::uu], c1[dIR3::uv] + c2[dIR3::uv],
      c1[dIR3::uw] + c2[dIR3::uw], c1[dIR3::vu] + c2[dIR3::vu],
      c1[dIR3::vv] + c2[dIR3::vv], c1[dIR3::vw] + c2[dIR3::vw],
      c1[dIR3::wu] + c2[dIR3::wu], c1[dIR3::wv] + c2[dIR3::wv],
      c1[dIR3::ww] + c2[dIR3::ww]};
} 



void equilibrium_a_field::build_interpolator_array(
    std::vector<std::unique_ptr<interpolator1d>>& interpolator_array,
    const narray_type& samples_array, const interpolator1d_factory* ifactory, size_t size_harmonics) {
  dblock_adapter sgrid(parser_->sgrid());
  
  std::transform(
      index_.begin(), index_.end(), interpolator_array.begin(), [&](size_t i) {
        std::slice mask_i(i, sgrid.size(), size_harmonics);
        narray_type data = samples_array[mask_i];
        return std::move( std::unique_ptr<interpolator1d>( ifactory->interpolate_data( sgrid, dblock_adapter(data) ) ) );
      });
}



//! Cached evaluation of trigonometric coefficients for internal Fourier series.
const equilibrium_a_field::cis_container_t& equilibrium_a_field::cached_cis_new(
    double theta, double zeta) const {
  thread_local double cached_theta = -1e6, cached_zeta = -1e6;
  thread_local cis_container_t cis_mn(mc_.size());
  if (theta != cached_theta || zeta != cached_zeta) {
    std::transform(
        index_.begin(), index_.end(), cis_mn.begin(),
        [&](size_t i) -> cis_container_t::value_type {
          double angle_mn = mc_[i] * theta - nc_[i] * zeta;
          return {std::cos(angle_mn), std::sin(angle_mn)};
        });
    cached_theta = theta;
    cached_zeta = zeta;
  }
  return cis_mn;
}



void equilibrium_a_field::build_harmonics(narray_type& new_m, narray_type& new_n, 
                                          const narray_type parser_m, const narray_type parser_n,
                                          narray_type& c_theta_val, narray_type& c_zeta_val) {
  //initialization of the auxiliars
  size_t size = parser_m.size();
  std::valarray<int> m_aux(2 * size * size);
  std::valarray<int> n_aux(2 * size * size);

  std::valarray<std::valarray<double>> C_theta_aux, C_zeta_aux;
  C_theta_aux.resize(2 * size * size, std::valarray<double>(parser_->ns()));
  C_zeta_aux.resize(2 * size * size, std::valarray<double>(parser_->ns()));

  std::vector<int> m_vec, n_vec;
  std::vector<std::valarray<double>> C_theta_vec, C_zeta_vec;


  // Fill the valarray's
  size_t k = 0;
  double norm_factor = (normalised_ ? this->m_factor() : 1);
  for (size_t i = 0; i<size; i++) {
      for (size_t j = 0; j<size; j++, k++) {
          m_aux[k] = parser_m[i] + parser_m[j];
          n_aux[k] = parser_n[i] + parser_n[j];

          for(size_t s=0; s<parser_->ns(); s++) {
              C_theta_aux[k][s] = -parser_->bsupumnc()[s*parser_->mnmax_nyq()+i] / norm_factor * parser_->gmnc()[s*parser_->mnmax_nyq()+j];
              C_zeta_aux[k][s] = -parser_->bsupvmnc()[s*parser_->mnmax_nyq()+i] / norm_factor * parser_->gmnc()[s*parser_->mnmax_nyq()+j];
          }

          k++;

          m_aux[k] = parser_m[i] - parser_m[j];
          n_aux[k] = parser_n[i] - parser_n[j];
          
          for(size_t s=0; s<parser_->ns(); s++) {
              C_theta_aux[k][s] = -parser_->bsupumnc()[s*parser_->mnmax_nyq()+i] / norm_factor * parser_->gmnc()[s*parser_->mnmax_nyq()+j];
              C_zeta_aux[k][s] = -parser_->bsupvmnc()[s*parser_->mnmax_nyq()+i] / norm_factor * parser_->gmnc()[s*parser_->mnmax_nyq()+j];
          }
      }
  }


  // Fill the vectors without rep
  m_vec.push_back(m_aux[0]);
  n_vec.push_back(n_aux[0]);
  C_theta_vec.push_back(C_theta_aux[0]);
  C_zeta_vec.push_back(C_zeta_aux[0]);
  bool flag;
  for (size_t i = 1; i<2*size*size; i++) {
      flag=1;
      for (size_t j = 0; j<m_vec.size(); j++) {
          if(m_vec[j]==m_aux[i] && n_vec[j]==n_aux[i] ) {
              C_theta_vec[j] += C_theta_aux[i];
              C_zeta_vec[j] += C_zeta_aux[i];
              flag=0;
          }
      }
      if (flag==1) {
          m_vec.push_back(m_aux[i]);
          n_vec.push_back(n_aux[i]);
          C_theta_vec.push_back(C_theta_aux[i]);
          C_zeta_vec.push_back(C_zeta_aux[i]);
      }
  }

  new_m.resize(m_vec.size());
  new_n.resize(m_vec.size());
  c_theta_val.resize(m_vec.size()*parser_->ns());
  c_zeta_val.resize(m_vec.size()*parser_->ns());
  for( size_t i = 0; i<m_vec.size(); i++ ) {
      new_m[i] = m_vec[i];
      new_n[i] = n_vec[i];

      for( size_t s=0; s<parser_->ns(); s++ ) {
          c_theta_val[s*m_vec.size()+i] = C_theta_vec[i][s];
          c_zeta_val[s*m_vec.size()+i] = C_zeta_vec[i][s];
      }
  }
}



void equilibrium_a_field::build_integrals(narray_type& ctheta_integ, narray_type& czeta_integ) {
  int index_00;
  double ds = 1.0/parser_->ns();
  

  for(size_t i=0; i<mc_.size(); i++) {
    if ( mc_[i]==0 && nc_[i]==0 ) {
      index_00 = i;
    }
  }
  

  ctheta_integ.resize(parser_->ns());
  czeta_integ.resize(parser_->ns());

  int i = 0;
  double ctheta_integ_aux = 0;
  double czeta_integ_aux = 0;
  
  for(double s=ds; s<=1; s=s+ds, i++) {
    ctheta_integ_aux = ctheta_integ_aux + ( (*ctheta_mn_[index_00])(s) + (*ctheta_mn_[index_00])(s-ds) ) * ds / 2;
    czeta_integ_aux = czeta_integ_aux + ( (*czeta_mn_[index_00])(s) + (*czeta_mn_[index_00])(s-ds) ) * ds / 2;

    ctheta_integ[i] = ctheta_integ_aux;
    czeta_integ[i] = czeta_integ_aux;
  }
}


}  // end namespace gyronimo.
