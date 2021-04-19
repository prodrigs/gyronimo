// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @parser_stellna.hh

#ifndef GYRONIMO_PARSER_STELLNA
#define GYRONIMO_PARSER_STELLNA

#include <string>
#include <iostream>
#include <valarray>

namespace gyronimo {

//! Parsing object for stellarator near-axis coordinates.
class parser_stellna {
 public:
  typedef std::valarray<double> narray_type;

  parser_stellna(const std::string& filename);
  ~parser_stellna() {};

  int Niota() const {return Niota_;};
  int field_periods() const {return field_periods_;};
  int n_phi() const {return n_phi_;};
  double R0() const {return R0_;};
  double eta_bar() const {return eta_bar_;};
  double iota() const {return iota_;};
  double axis_length() const {return axis_length_;};
  size_t axis_coeff_size() const {return axis_coeff_size_;};
  const narray_type& sigma() const {return sigma_;};
  const narray_type& Rcoeff() const {return Rcoeff_;};
  const narray_type& Zcoeff() const {return Zcoeff_;};
  const narray_type& normal() const {return normal_;};
  const narray_type& dldphi() const {return dldphi_;};
  const narray_type& torsion() const {return torsion_;};
  const narray_type& tangent() const {return tangent_;};
  const narray_type& binormal() const {return binormal_;};
  const narray_type& phi_grid() const {return phi_grid_;};
  const narray_type& curvature() const {return curvature_;};

 private:
  int Niota_, field_periods_, n_phi_;
  double R0_, eta_bar_, iota_, axis_length_;
  size_t axis_coeff_size_;
  narray_type Rcoeff_, Zcoeff_;
  narray_type phi_grid_, sigma_, curvature_, torsion_, dldphi_;
  narray_type tangent_, normal_, binormal_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_PARSER_HELENA
