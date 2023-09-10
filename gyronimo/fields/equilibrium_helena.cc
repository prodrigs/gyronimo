// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @equilibrium_helena.cc, this file is part of ::gyronimo::

#include <gyronimo/core/dblock.hh>
#include <gyronimo/fields/equilibrium_helena.hh>

namespace gyronimo{

equilibrium_helena::equilibrium_helena(
    const metric_helena *g, const interpolator2d_factory *ifactory)
    : IR3field_c1(std::abs(g->parser()->bmag()), 1.0, g),
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
  double chi = this->metric_->parser()->reduce_chi(position[IR3::v]);
  return {0.0, (*Bchi_)(s, chi), (*Bphi_)(s, chi)};
}
dIR3 equilibrium_helena::del_contravariant(
    const IR3& position, double time) const {
  double s = position[IR3::u];
  double chi = this->metric_->parser()->reduce_chi(position[IR3::v]);
  return {
      0.0, 0.0, 0.0, 
      Bchi_->partial_u(s, chi), Bchi_->partial_v(s, chi) , 0.0, 
      Bphi_->partial_u(s, chi), Bphi_->partial_v(s, chi) , 0.0};
}

}// end namespace gyronimo.
