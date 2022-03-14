// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @msphere_luhmann.cc

#include <cmath>
#include <gyronimo/fields/msphere_luhmann.hh>

namespace gyronimo {

msphere_luhmann::msphere_luhmann(
    double smooth_factor,
    double dipole_factor, double csheet_factor, double m_factor)
    : IR3field_c1(m_factor, 1.0, new metric_spherical(earth_radius)),
      c_bar_(0.001*csheet_factor/(earth_radius*m_factor)),
      d_bar_(dipole_factor/(earth_radius*m_factor)),
      idelta_(1.0/smooth_factor),
      metric_(nullptr) {
  metric_ = dynamic_cast<const metric_spherical*>(IR3field::metric());
}
msphere_luhmann::~msphere_luhmann() {
  if (metric_) delete metric_;
}
IR3 msphere_luhmann::contravariant(const IR3& position, double time) const {
  double r = position[IR3::u], r3 = r*r*r, r4 = r3*r;
  double cosv = std::cos(position[IR3::v]), sinv = std::sin(position[IR3::v]);
  double cosw = std::cos(position[IR3::w]), sinw = std::sin(position[IR3::w]);
  double tanh_factor = c_bar_*std::tanh(idelta_*r*cosv);
  double Bu = -2.0*d_bar_*cosv/r3 + tanh_factor*sinv*cosw;
  double Bv = -d_bar_*sinv/r4 + tanh_factor*cosv*cosw/r;
  double Bw = -tanh_factor*sinw/(r*sinv);
  return {Bu, Bv, Bw};
}
dIR3 msphere_luhmann::del_contravariant(
    const IR3& position, double time) const {
  double r = position[IR3::u], r2 = r*r, r3 = r2*r, r4 = r3*r, r5 = r4*r;
  double cosv = std::cos(position[IR3::v]), sinv = std::sin(position[IR3::v]);
  double cosw = std::cos(position[IR3::w]), sinw = std::sin(position[IR3::w]);
  double tanh_factor = c_bar_*std::tanh(idelta_*r*cosv);
  double sech_square = std::pow(std::cosh(idelta_*r*cosv), -2);
  double dBuu = cosv*(6*d_bar_/r4 + c_bar_*cosw*idelta_*sinv*sech_square);
  double dBuv = 2*d_bar_*sinv/r3 +
      cosw*(cosv*tanh_factor - c_bar_*idelta_*r*sinv*sinv*sech_square);
  double dBuw = -sinv*sinw*tanh_factor;
  double dBvu = 4*d_bar_*sinv/r5 + cosv*cosw*(
      c_bar_*cosv*idelta_*r*sech_square - tanh_factor)/r2;
  double dBvv = -cosv*d_bar_/r4 - cosw*sinv*(
      c_bar_*cosv*idelta_*r*sech_square + tanh_factor)/r;
  double dBvw = -cosv*sinw*tanh_factor/r;
  double dBwu = sinw*(
      tanh_factor/(r2*sinv) - c_bar_*cosv*idelta_*sech_square/(r*sinv));
  double dBwv = sinw*(
      c_bar_*idelta_*sech_square + cosv*tanh_factor/(r*sinv*sinv));
  double dBww = -cosw*tanh_factor/(r*sinv);
  return {dBuu, dBuv,  dBuw, dBvu, dBvv,  dBvw, dBwu, dBwv,  dBww};
}

} // end namespace gyronimo.
