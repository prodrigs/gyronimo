// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @parser_stellna.cc

#include <fstream>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/io_streaming.hh>
#include <gyronimo/parsers/parser_stellna.hh>

namespace gyronimo {

//! Reads and parses a STELLNA mapping file by name.
parser_stellna::parser_stellna(const std::string& filename) {
  std::ifstream input_stream;
  input_stream.open(filename);
  if (input_stream.rdstate() != std::ios_base::goodbit)
    error(__func__, __FILE__, __LINE__, "cannot open input file.", 1);

  input_stream >> R0_;
  input_stream >> axis_coeff_size_;
  Rcoeff_.resize(axis_coeff_size_);
  Zcoeff_.resize(axis_coeff_size_);
  input_stream >> Rcoeff_ >> Zcoeff_;
  input_stream >> eta_bar_;
  input_stream >> iota_;
  input_stream >> Niota_;
  input_stream >> axis_length_;
  input_stream >> field_periods_ >> n_phi_;
  for(auto p : {&phi_grid_, &sigma_,
      &curvature_, &torsion_, &dldphi_, &tangent_, &normal_, &binormal_}) {
    p->resize(n_phi_);
    input_stream >> (*p);
  }
}

} // end namespace gyronimo.
