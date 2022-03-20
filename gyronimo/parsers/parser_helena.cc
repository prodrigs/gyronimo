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

// @parser_helena.cc, this file is part of ::gyronimo::

#include <numbers>
#include <fstream>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/io_streaming.hh>
#include <gyronimo/parsers/parser_helena.hh>

namespace gyronimo {

//! Reads and parses a HELENA mapping file by name.
parser_helena::parser_helena(const std::string& filename) {
  std::ifstream input_stream;
  input_stream.open(filename);
  if (input_stream.rdstate() != std::ios_base::goodbit)
    error(__func__, __FILE__, __LINE__, "cannot open input file.", 1);

  input_stream >> npsi_; npsi_++;  // `HELENA` does't count the axis (s=0)...
  s_.resize(npsi_); input_stream >> s_;  // ... but stores it (shame on you)!
  q_.resize(npsi_); input_stream >> q_;

// Reads the dqs samples, which are stored in a rather peculiar form as:
// dqs[0], dqec, dq[1], ..., dq[npsi_ - 1]
  dqs_.resize(npsi_);
  input_stream >> dqec_;
  input_stream >> dqs_;
  std::swap(dqec_, dqs_[0]);

  curj_.resize(npsi_); input_stream >> curj_;
  input_stream >> dj0_ >> dje_;
  input_stream >> nchi_;
  is_symmetric_ = (nchi_ % 2 ? true : false);

// The angle \pi is stored by `HELENA` for symmetric equilibria, but for
// asymmetric ones 2*\pi is **not** stored. Weird, but it is the way it is...
  if (!is_symmetric_) nchi_++;
  chi_.resize(nchi_);

// Reads in the poloidal angles and adds the last one according to the symmetry.
  auto get_from_stream =
    [&input_stream]() {double x; input_stream >> x; return x;};
  std::generate(std::begin(chi_), (std::end(chi_) - 1), get_from_stream);
  if (is_symmetric_)  // symmetry exception...
    input_stream >> chi_[nchi_ - 1];
  else
    chi_[nchi_ - 1] = 2*std::numbers::pi;

// 2D fields are stored by `HELENA` without the row corresponding to the
// magnetic axis and without the column corresponding to the angle 2*\pi if
// not symmetric. Therefore, a special layout procedure is needed.
  this->layout_2d_field(input_stream, gmh11_);
  this->layout_2d_field(input_stream, gmh12_);

  input_stream >> cpsurf_ >> radius_;
  this->layout_2d_field(input_stream, gmh33_);

  input_stream >> raxis_;  // sadly, HELENA always prints 1.0 here.
  p0_.resize(npsi_); input_stream >> p0_;
  input_stream >> dp0_ >> dpe_;
  rbphi_.resize(npsi_); input_stream >> rbphi_;
  input_stream >> drbphi0_ >> drbphie_;

// Reads in vaccum information and add the last angular value.
  vx_.resize(nchi_);
  std::generate(std::begin(vx_), (std::end(vx_) - 1), get_from_stream);
  if (is_symmetric_)  // symmetry exception...
    input_stream >> vx_[nchi_ - 1];
  else
    vx_[nchi_ - 1] = vx_[0];
  vy_.resize(nchi_);
  std::generate(std::begin(vy_), (std::end(vy_) - 1), get_from_stream);
  if (is_symmetric_)  // symmetry exception...
    input_stream >> vy_[nchi_ - 1];
  else
    vy_[nchi_ - 1] = vy_[0];

  input_stream >> eps_;

  this->layout_2d_field(input_stream, x_);
  this->layout_2d_field(input_stream, y_);

  input_stream >> rmag_ >> bmag_;
  input_stream.close();

  rgeo_ = radius_/eps_*rmag_;
  this->build_auxiliar_data();
}
 
//! Builds up some useful flux-function 2D arrays.
/*
    Each flux-function value is copied over all poloidal samples along the flux
    surface (or row) in order to allow arithmetic operations with other non
    flux-function 2D arrays.
*/
void parser_helena::build_auxiliar_data() {
  f_.resize(npsi_*nchi_);
  F_.resize(npsi_*nchi_);
  qoF_.resize(npsi_*nchi_);

// Replicates 1d radial profiles over the 2d grid:
  for (size_t row = 0; row < npsi_; row++) {
    std::slice flux_surface(row*nchi_, nchi_, 1);
    F_[flux_surface] = rbphi_[row];
    f_[flux_surface] = 2.0*cpsurf_*s_[row];
    qoF_[flux_surface] = q_[row]/rbphi_[row];
  }

  J_ = f_*qoF_*gmh33_;
  covariant_g33_ = gmh33_;
  covariant_g22_ = qoF_*qoF_*gmh33_*gmh11_;
  covariant_g12_ = -qoF_*qoF_*f_*gmh33_*gmh12_;
  covariant_g11_ = (1.0 + qoF_*qoF_*gmh12_*gmh12_*gmh33_)*f_*f_/gmh11_;

// Corrects the indeterminate ratio at the axis:
  covariant_g11_[std::slice(0, nchi_, 1)] =
    this->axis_extrapolation(covariant_g11_);

  covariant_B1_ = -f_*qoF_*gmh12_;
  covariant_B2_ = qoF_*gmh11_;
  covariant_B3_ = F_;

  contravariant_B1_ = narray_type(0.0, nchi_*npsi_);
  contravariant_B2_ = 1.0/(qoF_*gmh33_);  // Note: f_/J_ is indeterminate at 0!
  contravariant_B3_ = F_/gmh33_;
}

//! Extrapolates the axis' row value from neighbouring points.
double parser_helena::axis_extrapolation(const narray_type& array) {
  double a1 = array[nchi_], a2 = array[2*nchi_], a3 = array[3*nchi_];
  double ds0 = s_[1] - s_[0], ds2 = s_[2] - s_[1], ds3 = s_[3] - s_[1];
  return (a1*(ds0 + ds2)*(ds2 - ds3)*(ds0 + ds3) + 
    ds0*(-(a3*ds2*(ds0 + ds2)) + a2*ds3*(ds0 + ds3)))/(ds2*(ds2 - ds3)*ds3);
}

//! Adds the axis row and the final angle column (if asym) to a 2D field.
void parser_helena::layout_2d_field(
    std::istream& input_stream, narray_type& composed_array) {
  std::slice axis_row(0, nchi_, 1);
  size_t stored_nchi = (is_symmetric_ ? nchi_ : nchi_ - 1);
  std::gslice inner_submatrix(nchi_, {npsi_ - 1, stored_nchi}, {nchi_, 1});

  narray_type array_from_hmap(stored_nchi*(npsi_ - 1));
  input_stream >> array_from_hmap;

  composed_array.resize(npsi_*nchi_);
  composed_array[inner_submatrix] = array_from_hmap;
  composed_array[axis_row] = this->axis_extrapolation(composed_array);

  if (!is_symmetric_) {  // symmetry exception...
    std::slice first_column(0, npsi_, nchi_);
    std::slice last_column(nchi_ - 1, npsi_, nchi_);
    composed_array[last_column] = composed_array[first_column];
  }
}

} // end namespace gyronimo.
