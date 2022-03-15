// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @field_line.hh

#ifndef GYRONIMO_FIELD_LINE
#define GYRONIMO_FIELD_LINE

#include <array>
#include <gyronimo/fields/IR3field.hh>

namespace gyronimo {

//! Field-line equations of motion parametrised by length over `Lref`.
class field_line {
 public:
  typedef std::array<double, 3> state;
  field_line(
      const IR3field* field, double Lref) : field_(field), Lref_(Lref) {};
  ~field_line() {};
  state operator()(const state& x, double t) const {
      return IR3(Lref_*field_->contravariant_versor(x, t));};
  const IR3field* field() const {return field_;};
  double Lref() const {return Lref_;};
 private:
  const IR3field* field_;
  const double Lref_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_FIELD_LINE
