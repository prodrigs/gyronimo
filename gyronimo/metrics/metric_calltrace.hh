// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2023 Paulo Rodrigues.

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

// @metric_calltrace.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_CALLTRACE
#define GYRONIMO_METRIC_CALLTRACE

#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! Calltrace layer over objects derived from `metric_covariant`.
/*!
    Defines a decorator layer to trace calls to each member function declared in
    the `metric_covariant` interface. Because it derives from its template
    argument, it can be used to replace pointers or references to objects of
    that type in client code. The cache object is constructed with the same
    constructor call as its template argument.
*/
template<typename T> requires std::derived_from<T, metric_covariant>
class metric_calltrace : public T {
 public:
  template<typename... Args> metric_calltrace(Args... args) : T(args...) {};
  virtual ~metric_calltrace() {};
  virtual SM3 operator()(const IR3& q) const override;
  virtual dSM3 del(const IR3& q) const override;
  virtual double jacobian(const IR3& q) const override;
  virtual IR3 del_jacobian(const IR3& q) const override;
  virtual SM3 inverse(const IR3& q) const override;
  virtual dSM3 del_inverse(const IR3& q) const override;
  virtual ddIR3 christoffel_first_kind(const IR3& q) const override;
  virtual ddIR3 christoffel_second_kind(const IR3& q) const override;
  virtual IR3 inertial_force(const IR3& q, const IR3& dot_q) const override;
  virtual IR3 to_covariant(const IR3& B, const IR3& q) const override;
  virtual IR3 to_contravariant(const IR3& B, const IR3& q) const override;
};

template<typename T> requires std::derived_from<T, metric_covariant>
SM3 metric_calltrace<T>::operator()(const IR3& q) const {
  std::cout << "metric_covariant::operator()(q)\n";
  return T::operator()(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
dSM3 metric_calltrace<T>::del(const IR3& q) const {
  std::cout << "metric_covariant::del(q)\n";
  return T::del(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
double metric_calltrace<T>::jacobian(const IR3& q) const {
  std::cout << "metric_covariant::jacobian(q)\n";
  return T::jacobian(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_calltrace<T>::del_jacobian(const IR3& q) const {
  std::cout << "metric_covariant::del_jacobian(q)\n";
  return T::del_jacobian(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
SM3 metric_calltrace<T>::inverse(const IR3& q) const {
  std::cout << "metric_covariant::inverse(q)\n";
  return T::inverse(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
dSM3 metric_calltrace<T>::del_inverse(const IR3& q) const {
  std::cout << "metric_covariant::del_inverse(q)\n";
  return T::del_inverse(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
ddIR3 metric_calltrace<T>::christoffel_first_kind(const IR3& q) const {
  std::cout << "metric_covariant::christoffel_first_kind(q)\n";
  return T::christoffel_first_kind(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
ddIR3 metric_calltrace<T>::christoffel_second_kind(const IR3& q) const {
  std::cout << "metric_covariant::christoffel_second_kind(q)\n";
  return T::christoffel_second_kind(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_calltrace<T>::inertial_force(const IR3& q, const IR3& dot_q) const {
  std::cout << "metric_covariant::inertial_force(q, dot_q)\n";
  return T::inertial_force(q, dot_q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_calltrace<T>::to_covariant(const IR3& B, const IR3& q) const {
  std::cout << "metric_covariant::to_covariant(B, q)\n";
  return T::to_covariant(B, q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_calltrace<T>::to_contravariant(const IR3& B, const IR3& q) const {
  std::cout << "metric_covariant::to_contravariant(B, q)\n";
  return T::to_contravariant(B, q);
}

}  // namespace gyronimo

#endif  // GYRONIMO_METRIC_CALLTRACE
