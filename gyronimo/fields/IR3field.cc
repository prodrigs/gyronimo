// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @IR3field.cc

#include <cmath>
#include <limits>
#include <gyronimo/core/error.hh>
#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/core/contraction.hh>

namespace gyronimo {

IR3field::IR3field(double m_factor, double t_factor, const metric_covariant* g)
    : m_factor_(m_factor), t_factor_(t_factor), covariant_metric_(g) {
  if (!g) error(__func__, __FILE__, __LINE__, "invalid metric pointer.", 1);
  if (t_factor < std::numeric_limits<double>::epsilon())
    error(__func__, __FILE__, __LINE__, "non_positive t_factor.", 1);
  if (m_factor < std::numeric_limits<double>::epsilon())
    error(__func__, __FILE__, __LINE__, "non-positive m_factor.", 1);
}
IR3 IR3field::covariant(const IR3& position, double time) const {
  IR3 A = this->contravariant(position, time);
  return covariant_metric_->to_covariant(A, position);
}
double IR3field::magnitude(const IR3& position, double time) const {
  IR3 A = this->contravariant(position, time);
  IR3 B = covariant_metric_->to_covariant(A, position);
  return std::sqrt(inner_product(A, B));
}
IR3 IR3field::covariant_versor(const IR3& position, double time) const {
  double imagnitude = 1.0/this->magnitude(position, time);
  IR3 A = this->covariant(position, time);
  return {imagnitude*A[IR3::u], imagnitude*A[IR3::v], imagnitude*A[IR3::w]};
}
IR3 IR3field::contravariant_versor(const IR3& position, double time) const {
  double imagnitude = 1.0/this->magnitude(position, time);
  IR3 A = this->contravariant(position, time);
  return {imagnitude*A[IR3::u], imagnitude*A[IR3::v], imagnitude*A[IR3::w]};
}

} // end namespace gyronimo.
