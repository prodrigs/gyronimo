// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues.

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

// @parser_castor.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_PARSER_CASTOR
#define GYRONIMO_PARSER_CASTOR

#include <string>
#include <valarray>

namespace gyronimo {

//! Parsing object for `CASTOR` ceig output files.
/*!
    Reads and parses a ceig output file produced by the MHD eigenvalue code
    `CASTOR` [K. Kerner *et al*., J. Comput. Phys. **142**, 271 (1998)]. The
    variables v1, v2, v3, a2, a3, t, and rho are provided as defined in the
    published paper. As an exception, a1 is provided as the true first covariant
    component of the potential vector, not the product -i*a1 stored in the ceig
    file.
*/

class parser_castor {
 public:
  typedef std::valarray<double> narray_type;

  parser_castor(const std::string &filename);
  ~parser_castor() {};

  size_t n_psi() const {return n_psi_;};
  size_t n_harm() const {return n_harm_;};
  double n_tor() const {return n_tor_;};
  double eigenvalue_real() const {return eigenvalue_real_;};
  double eigenvalue_imag() const {return eigenvalue_imag_;};
  const narray_type& s() const {return s_;};
  const narray_type& m() const {return m_;};
  const narray_type& t_real() const {return t_real_;};
  const narray_type& v1_real() const {return v1_real_;};
  const narray_type& v2_real() const {return v2_real_;};
  const narray_type& v3_real() const {return v3_real_;};
  const narray_type& a1_real() const {return a1_real_;};
  const narray_type& a2_real() const {return a2_real_;};
  const narray_type& a3_real() const {return a3_real_;};
  const narray_type& rho_real() const {return rho_real_;};
  const narray_type& t_imag() const {return t_imag_;};
  const narray_type& v1_imag() const {return v1_imag_;};
  const narray_type& v2_imag() const {return v2_imag_;};
  const narray_type& v3_imag() const {return v3_imag_;};
  const narray_type& a1_imag() const {return a1_imag_;};
  const narray_type& a2_imag() const {return a2_imag_;};
  const narray_type& a3_imag() const {return a3_imag_;};
  const narray_type& rho_imag() const {return rho_imag_;};

 private:
  double n_tor_;
  size_t n_psi_, n_harm_;
  double eigenvalue_real_, eigenvalue_imag_;
  narray_type s_, m_;
  narray_type rho_real_, rho_imag_, t_real_, t_imag_;
  narray_type v1_real_, v1_imag_, v2_real_, v2_imag_, v3_real_, v3_imag_;
  narray_type a1_real_, a1_imag_, a2_real_, a2_imag_, a3_real_, a3_imag_;

  void initialise_variable_chunk(
      std::ifstream& input_stream, narray_type& real, narray_type& imag);
};

} // end namespace gyronimo.

#endif // GYRONIMO_PARSER_CASTOR
