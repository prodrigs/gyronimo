// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @parser_helena.hh

#ifndef GYRONIMO_PARSER_HELENA
#define GYRONIMO_PARSER_HELENA

#include <string>
#include <iostream>
#include <valarray>

namespace gyronimo {

//! Parsing object for `HELENA` mapping files.
/*!
    Reads and parses a mapping file produced by `HELENA`, a Grad-Shafanov
    equilibrium code [G. Huysmans *et al*., Int. J.  Mod. Phys. C 2 **371**,
    (1991)]. The mapping file is usually produced as text file named `fort.12`.
    Parsed data is accessed via public member functions, some of them returning
    `parser_helena::narray_type` arrays. Useful data is computed from the parsed
    one and made available by the same mechanism: Contravariant metric-tensor
    components are normalised to @f$R_0^{-2}@f$, covariant to @f$R_0^2@f$,
    whilst contravariant and covariant magnetic field components are normalised,
    respectively, to @f$B_0/R_0@f$ and @f$B_0 R_0@f$, where @f$R_0@f$ and
    @f$B_0@f$ are the location of the magnetic axis and the field magnitude
    there [rmag() and bmag(), both in SI units]. The coordinate set is @f$\{s,
    \chi, \phi \}@f$, where @f$0 \leq s \leq 1@f$ is the square root of the
    poloidal flux normalised to its boundary value (in Wb/rad), @f$\chi@f$ is a
    counterclockwise angle such that @f$B^\phi=q B^\chi@f$, with @f$\chi=0@f$
    the low-field side midplane, @f$\phi@f$ is the clockwise toroidal angle as
    seen from above the torus, and @f$q@f$ the safety factor.
*/

class parser_helena {
 public:
  typedef std::valarray<double> narray_type;

  parser_helena(const std::string& filename);
  ~parser_helena() {};

  bool is_symmetric() const {return is_symmetric_;};
  size_t npsi() const {return npsi_;};  // Radial-grid size.
  size_t nchi() const {return nchi_;};  // Poloidal-grid size.
  double cpsurf() const {return cpsurf_;};  // Edge flux over @f$B_0R_0^2@f$.
  double radius() const {return radius_;};  // Ratio @f$a/R_0@f$.
  double eps() const {return eps_;};  // Ratio @f$a/R_{geo}@f$.
  double rgeo() const {return rgeo_;};  // Grid origin @f$R_{geo}@f$ (m).
  double rmag() const {return rmag_;};  // Axis position @f$R_0@f$ (m).
  double bmag() const {return bmag_;};  // Axis magnetic field @f$B_0@f$.
  double dqec() const {return dqec_;};
  double dj0() const {return dj0_;};
  double dje() const {return dje_;};
  double dp0() const {return dp0_;};
  double dpe() const {return dpe_;};
  double drbphi0() const {return drbphi0_;};
  double drbphie() const {return drbphie_;};
  const narray_type& s() const {return s_;};  // Radial-grid nodes.
  const narray_type& q() const {return q_;};  // Safety-factor samples.
  const narray_type& p0() const {return p0_;};  // @f$\mu_0 p/B_0^2@f$.
  const narray_type& dqs() const {return dqs_;};
  const narray_type& chi() const {return chi_;};  // Poloidal-grid samples.
  const narray_type& curj() const {return curj_;};  // @f$<J>/J_0@f$.
  const narray_type& rbphi() const {return rbphi_;}; // @f$B_\phi/(R_0 B_0)@f$.
  const narray_type& gmh11() const {return gmh11_;};
  const narray_type& gmh12() const {return gmh12_;};
  const narray_type& gmh33() const {return gmh33_;};
  const narray_type& vx() const {return vx_;};
  const narray_type& vy() const {return vy_;};
  const narray_type& x() const {return x_;};  // @f$(R - R_{geo})/a@f$ at grid.
  const narray_type& y() const {return y_;};  // @f$(Z - Z_{geo})/a@f$ at grid.
  const narray_type& f() const {return f_;};
  const narray_type& F() const {return F_;};
  const narray_type& qoF() const {return qoF_;};
  const narray_type& covariant_g11() const {return covariant_g11_;};
  const narray_type& covariant_g12() const {return covariant_g12_;};
  const narray_type& covariant_g22() const {return covariant_g22_;};
  const narray_type& covariant_g33() const {return covariant_g33_;};
  const narray_type& covariant_B1() const {return covariant_B1_;};
  const narray_type& covariant_B2() const {return covariant_B2_;};
  const narray_type& covariant_B3() const {return covariant_B3_;};
  const narray_type& contravariant_B2() const {return contravariant_B2_;};
  const narray_type& contravariant_B3() const {return contravariant_B3_;};

 private:
  bool is_symmetric_;
  size_t npsi_, nchi_;
  double cpsurf_, radius_, raxis_, eps_, rgeo_, rmag_, bmag_;
  double dqec_, dj0_, dje_, dp0_, dpe_, drbphi0_, drbphie_;
  narray_type s_, q_, p0_, dqs_, chi_, curj_;
  narray_type rbphi_, gmh11_, gmh12_, gmh33_;
  narray_type vx_, vy_, x_, y_;
  narray_type f_, F_, qoF_, J_;
  narray_type covariant_g11_, covariant_g12_;
  narray_type covariant_g22_, covariant_g33_;
  narray_type covariant_B1_, covariant_B2_, covariant_B3_;
  narray_type contravariant_B1_, contravariant_B2_, contravariant_B3_;

  void build_auxiliar_data();
  double axis_extrapolation(const narray_type& array);
  void layout_2d_field(std::istream& input_stream, narray_type& composed_array);
};

} // end namespace gyronimo.

#endif // GYRONIMO_PARSER_HELENA
