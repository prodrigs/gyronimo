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

// @IR3field_c1_calltrace.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_IR3FIELD_C1_CALLTRACE
#define GYRONIMO_IR3FIELD_C1_CALLTRACE

#include <gyronimo/fields/IR3field_c1.hh>

namespace gyronimo {

//! Cache layer over objects derived from `IR3field_c1`.
/*!
    Defines a decorator layer to trace calls to each member function declared in
    the `IR3field_c1` interface. Because it derives from its template argument,
    it can be used to replace pointers or references to objects of that type in
    client code. The cache object is constructed with the same constructor call
    as its template argument.
*/
template<typename T> requires std::derived_from<T, IR3field_c1>
class IR3field_c1_calltrace : public T {
 public:
  template<typename... Args>
  IR3field_c1_calltrace(Args... args) : T(args...) {};
  virtual ~IR3field_c1_calltrace() {};
  virtual IR3 contravariant(const IR3& q, double t) const override;
  virtual IR3 covariant(const IR3& q, double t) const override;
  virtual double magnitude(const IR3& q, double t) const override;
  virtual IR3 covariant_versor(const IR3& q, double t) const override;
  virtual IR3 contravariant_versor(const IR3& q, double t) const override;
  virtual dIR3 del_contravariant(const IR3& q, double t) const override;
  virtual IR3 partial_t_contravariant(const IR3& q, double t) const override;
  virtual IR3 del_magnitude(const IR3& q, double t) const override;
  virtual double partial_t_magnitude(const IR3& q, double t) const override;
  virtual dIR3 del_covariant(const IR3& q, double t) const override;
  virtual IR3 partial_t_covariant(const IR3& q, double t) const override;
  virtual IR3 curl(const IR3& q, double t) const override;
};

template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::contravariant(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::contravariant(q, t)\n";
  return T::contravariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::covariant(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::covariant(q, t)\n";
  return T::covariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
double IR3field_c1_calltrace<T>::magnitude(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::magnitude(q, t)\n";
  return T::magnitude(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::covariant_versor(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::covariant_versor(q, t)\n";
  return T::covariant_versor(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::contravariant_versor(
    const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::contravariant_versor(q, t)\n";
  return T::contravariant_versor(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
dIR3 IR3field_c1_calltrace<T>::del_contravariant(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::del_contravariant(q, t)\n";
  return T::del_contravariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::partial_t_contravariant(
    const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::partial_t_contravariant(q, t)\n";
  return T::partial_t_contravariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::del_magnitude(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::del_magnitude(q, t)\n";
  return T::del_magnitude(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
double IR3field_c1_calltrace<T>::partial_t_magnitude(
    const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::partial_t_magnitude(q, t)\n";
  return T::partial_t_magnitude(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
dIR3 IR3field_c1_calltrace<T>::del_covariant(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::del_covariant(q, t)\n";
  return T::del_covariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::partial_t_covariant(
    const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::partial_t_covariant(q, t)\n";
  return T::partial_t_covariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_calltrace<T>::curl(const IR3& q, double t) const {
  std::cout << "IR3fiedl_c1::curl(q, t)\n";
  return T::curl(q, t);
}

}  // namespace gyronimo

#endif  // GYRONIMO_IR3FIELD_C1_CALLTRACE
