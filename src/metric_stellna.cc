// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @metric_stellna.cc

#include <numbers>
#include <gyronimo/metrics/metric_stellna.hh>

namespace gyronimo{

metric_stellna::metric_stellna(
    const parser_stellna *parser, const interpolator1d_factory *ifactory)
    : parser_(parser),
      sigma_(nullptr), curvature_(nullptr), torsion_(nullptr), dldphi_(nullptr),
      phi_modulus_factor_(2*std::numbers::pi/parser->field_periods()) {
  dblock_adapter phi_grid(parser->phi_grid());
  sigma_ = ifactory->interpolate_data(
      phi_grid, dblock_adapter(parser->sigma()));
  dldphi_ = ifactory->interpolate_data(
      phi_grid, dblock_adapter(parser->dldphi()));
  torsion_ = ifactory->interpolate_data(
      phi_grid, dblock_adapter(parser->torsion()));
  curvature_ = ifactory->interpolate_data(
      phi_grid, dblock_adapter(parser->curvature()));
}
metric_stellna::~metric_stellna() {
  if(curvature_) delete curvature_;
  if(torsion_) delete torsion_;
  if(dldphi_) delete dldphi_;
  if(sigma_) delete sigma_;
}
SM3 metric_stellna::operator()(const IR3& position) const {
  double phi = std::fmod(position[IR3::w], phi_modulus_factor_);
  double r = position[IR3::u], theta = position[IR3::v];
  double coso = std::cos(theta), sino = std::sin(theta);
  double l_prime = (*dldphi_)(phi);
  double k = (*curvature_)(phi), k_prime = (*curvature_).derivative(phi);
  double sigma = (*sigma_)(phi), sigma_prime = (*sigma_).derivative(phi);
  double eta_over_k = parser_->eta_bar()/k;
  double eta_over_k_squared = eta_over_k*eta_over_k;
  double guu = std::pow(eta_over_k*coso, 2) +
          std::pow((sino + sigma*coso)/eta_over_k, 2);
  double guv = -r*eta_over_k_squared*sino*coso +
      r/eta_over_k_squared*(sino + sigma*coso)*(coso - sigma*sino);
  double guw = r*(sigma_prime*coso*(sino + sigma*coso) +
      k_prime/k*(
          std::pow(sino + sigma*coso, 2) -
          std::pow(eta_over_k_squared*coso, 2)))/eta_over_k_squared;
  double gvv = r*r*(std::pow(eta_over_k*sino, 2) +
          std::pow((coso - sigma*sino)/eta_over_k, 2));
  double factor_sigma_k = sigma*k_prime/k + 0.5*sigma_prime;
  double gvw = r*r*(l_prime*(*torsion_)(phi) + (
          0.5*sigma_prime +
          factor_sigma_k*(coso*coso - sino*sino) +
          coso*sino*(k_prime/k*(1.0 - sigma*sigma) - sigma*sigma_prime)
      )/eta_over_k_squared + coso*sino*eta_over_k_squared*k_prime/k);
  double gww = l_prime*l_prime*(1.0 - 2.0*parser_->eta_bar()*r*coso);
  return {guu, guv, guw, gvv, gvw, gww};
}
dSM3 metric_stellna::del(const IR3& position) const {
  double phi = std::fmod(position[IR3::w], phi_modulus_factor_);
  double r = position[IR3::u], theta = position[IR3::v];
  double coso = std::cos(theta), sino = std::sin(theta);
  double l_prime = (*dldphi_)(phi), l_prime_prime=(*dldphi_).derivative(phi);
  double k = (*curvature_)(phi), k_prime = (*curvature_).derivative(phi),
      k_prime_prime = (*curvature_).derivative2(phi);
  double sigma = (*sigma_)(phi), sigma_prime = (*sigma_).derivative(phi),
      sigma_prime_prime = (*sigma_).derivative2(phi);
  double tau = (*torsion_)(phi), tau_prime = (*torsion_).derivative(phi);
  double eta_over_k = parser_->eta_bar()/k;
  double eta_over_k_squared = eta_over_k*eta_over_k;
  double factor_a = 1.0 + eta_over_k_squared*eta_over_k_squared - sigma*sigma;
  double factor_b = sigma*sigma_prime - factor_a*k_prime/k;
  double cos2o = coso*coso - sino*sino;
  double d_u_guu = 0.0;
  double d_v_guu = 2.0*(sino + sigma*coso)*(
      coso - sigma*sino)/eta_over_k_squared - 2.0*eta_over_k_squared*coso*sino;
  double d_w_guu =(
      2.0*sino*sino*k_prime/k +
      coso*sino*(4.0*sigma*k_prime/k + 2.0*sigma_prime) +
      coso*coso*2.0*(sigma*sigma_prime - k_prime/k*(
          eta_over_k_squared*eta_over_k_squared - sigma*sigma))
      )/eta_over_k_squared;
  double d_u_guv = -coso*sino*eta_over_k_squared +
      (sino + sigma*coso)*(coso - sigma*sino)/eta_over_k_squared;
  double d_v_guv = r*((sino*sino - coso*coso)*(
      sigma*sigma + eta_over_k_squared*eta_over_k_squared - 1.0) -
          4.0*coso*sino*sigma)/eta_over_k_squared;
  double d_w_guv = r*(
      cos2o*(2.0*sigma*k_prime/k + sigma_prime) -
      2.0*coso*sino*factor_b)/eta_over_k_squared;
  double d_u_guw = (sino*sino*k_prime/k +
      coso*sino*(2.0*sigma*k_prime/k + sigma_prime) +
      coso*coso*(sigma*sigma_prime +
          k_prime/k*(sigma*sigma - eta_over_k_squared*eta_over_k_squared))
      )/eta_over_k_squared;
  double d_v_guw = r*(cos2o*(2.0*sigma*k_prime/k + sigma_prime) -
      2.0*coso*sino*factor_b)/eta_over_k_squared;
  double factor_k_prime2 = k_prime/k*k_prime/k + k_prime_prime/k;
  double d_w_guw = r*(
      sino*sino*factor_k_prime2 +
      coso*sino*(2.0*sigma*factor_k_prime2 + 4.0*k_prime*sigma_prime/k +
      sigma_prime_prime) + coso*coso*(
          k_prime/k*k_prime/k*(
              3.0*eta_over_k_squared*eta_over_k_squared + sigma*sigma) +
          4.0*sigma*sigma_prime*k_prime/k + sigma_prime*sigma_prime +
          sigma*sigma_prime_prime + k_prime_prime/k*(
              sigma*sigma - eta_over_k_squared*eta_over_k_squared))
      )/eta_over_k_squared;
  double d_u_gvv = 2.0*r*(coso*coso - 2.0*coso*sino*sigma +
      sino*sino*(eta_over_k_squared*eta_over_k_squared + sigma*sigma)
      )/eta_over_k_squared;
  double d_v_gvv = 2.0*r*r*(sigma*(sino*sino - coso*coso) +
      coso*sino*(sigma*sigma + eta_over_k_squared*eta_over_k_squared - 1.0)
      )/eta_over_k_squared;
  double d_w_gvv = 2.0*r*r*(coso*coso*k_prime/k -
      coso*sino*(2.0*sigma*k_prime/k+ sigma_prime) +
      sino*sino*(sigma*sigma_prime -
          (eta_over_k_squared*eta_over_k_squared - sigma*sigma)*k_prime/k)
      )/eta_over_k_squared;
  double d_u_gvw = 2.0*r*l_prime*(*torsion_)(phi) + r*(
      sigma_prime + cos2o*(2.0*sigma*k_prime/k + sigma_prime) -
      2.0*coso*sino*factor_b)/eta_over_k_squared;
  double d_v_gvw = -r*r*(cos2o*factor_b +
          2.0*coso*sino*(2.0*sigma*k_prime/k + sigma_prime))/eta_over_k_squared;
  double d_w_gvw = r*r*(eta_over_k_squared*(
          l_prime*tau_prime +
          tau*l_prime_prime) +
      0.5*sigma_prime_prime + sigma_prime*k_prime/k +
      0.5*cos2o*(
          2.0*sigma*(k_prime/k*k_prime/k + k_prime_prime/k) +
          4.0*sigma_prime*k_prime/k + sigma_prime_prime) -
      coso*sino*(
          k_prime/k*k_prime/k*(
              sigma*sigma + 3.0*eta_over_k_squared*eta_over_k_squared - 1.0) +
          4.0*sigma*sigma_prime*k_prime/k + sigma_prime*sigma_prime +
          sigma*sigma_prime_prime - k_prime_prime/k*factor_a)
      )/eta_over_k_squared;
  double d_u_gww = -2.0*parser_->eta_bar()*l_prime*l_prime*coso;
  double d_v_gww = l_prime*l_prime*2.0*parser_->eta_bar()*r*sino;
  double d_w_gww = 2.0*l_prime*l_prime_prime*(
      1.0 - 2.0*parser_->eta_bar()*r*coso);
  return {
      d_u_guu, d_v_guu, d_w_guu,
      d_u_guv, d_v_guv, d_w_guv,
      d_u_guw, d_v_guw, d_w_guw,
      d_u_gvv, d_v_gvv, d_w_gvv,
      d_u_gvw, d_v_gvw, d_w_gvw,
      d_u_gww, d_v_gww, d_w_gww};
}

} // end namespace gyronimo
