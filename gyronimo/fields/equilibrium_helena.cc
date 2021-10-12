// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @equilibrium_helena.cc

#include <gyronimo/core/dblock.hh>
#include <gyronimo/core/transpose.hh>
#include <gyronimo/fields/equilibrium_helena.hh>

namespace gyronimo{

equilibrium_helena::equilibrium_helena(
    const metric_helena *g, const interpolator2d_factory *ifactory)
    : IR3field_c1(g->parser()->bmag(), 1.0, g),
      metric_(g), Bchi_(nullptr), Bphi_(nullptr) {
  const parser_helena *p = metric_->parser();
  dblock_adapter s_range(p->s()), chi_range(p->chi());
  Bchi_ = ifactory->interpolate_data(
      s_range, chi_range,
      dblock_adapter(parser_helena::narray_type(p->contravariant_B2()/R0())));
  Bphi_ = ifactory->interpolate_data(
      s_range, chi_range,
      dblock_adapter(parser_helena::narray_type(p->contravariant_B3()/R0())));
}
equilibrium_helena::~equilibrium_helena() {
  if(Bchi_) delete Bchi_;
  if(Bphi_) delete Bphi_;
}
IR3 equilibrium_helena::contravariant(const IR3& position, double time) const {
  double s = position[IR3::u];
  double chi = this->metric_->reduce_chi(position[IR3::v]);
  return {0.0, (*Bchi_)(s, chi), (*Bphi_)(s, chi)};
}
dIR3 equilibrium_helena::del_contravariant(
    const IR3& position, double time) const {
  double s = position[IR3::u];
  double chi = this->metric_->reduce_chi(position[IR3::v]);
  return {
      0.0, 0.0, 0.0, 
	  Bchi_->partial_u(s, chi), Bchi_->partial_v(s, chi) , 0.0, 
	  Bphi_->partial_u(s, chi), Bphi_->partial_v(s, chi) , 0.0};
}

}// end namespace gyronimo.
