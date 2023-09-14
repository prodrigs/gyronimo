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

// @morphism_calltrace.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_CALLTRACE
#define GYRONIMO_MORPHISM_CALLTRACE

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Cache layer over objects derived from `morphism`.
/*!
    Defines a decorator layer to trace calls to each member function declared in
    the `morphism` interface. Because it derives from its template argument, it
    can be used to replace pointers or references to objects of that type in
    client code. The cache object is constructed with the same constructor call
    as its template argument.
*/
template<typename T> requires std::derived_from<T, morphism>
class morphism_calltrace : public T {
 public:
  template<typename... Args> morphism_calltrace(Args... args) : T(args...) {};
  virtual ~morphism_calltrace() {};
  virtual IR3 operator()(const IR3& q) const override;
  virtual IR3 inverse(const IR3& x) const override;
  virtual dIR3 del(const IR3& q) const override;
  virtual ddIR3 ddel(const IR3& q) const override;
  virtual double jacobian(const IR3& q) const override;
  virtual dIR3 del_inverse(const IR3& q) const override;
  virtual dIR3 tan_basis(const IR3& q) const override;
  virtual dIR3 dual_basis(const IR3& q) const override;
  virtual IR3 to_covariant(const IR3& A, const IR3& q) const override;
  virtual IR3 to_contravariant(const IR3& A, const IR3& q) const override;
  virtual IR3 from_covariant(const IR3& A, const IR3& q) const override;
  virtual IR3 from_contravariant(const IR3& A, const IR3& q) const override;
  virtual IR3 translation(const IR3& q, const IR3& delta) const override;
 private:
  static inline IR3 un_init_ = {123456789., 987654321., 192837465.};
};

template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_calltrace<T>::operator()(const IR3& q) const {
  std::cout << "morphism::operator()(q)\n";
  return T::operator()(q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_calltrace<T>::inverse(const IR3& x) const {
  std::cout << "morphism::inverse(x)\n";
  return T::inverse(x);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_calltrace<T>::del(const IR3& q) const {
  std::cout << "morphism::del(q)\n";
  return T::del(q);
}
template<typename T> requires std::derived_from<T, morphism>
ddIR3 morphism_calltrace<T>::ddel(const IR3& q) const {
  std::cout << "morphism::ddel(q)\n";
  return T::ddel(q);
}
template<typename T> requires std::derived_from<T, morphism>
double morphism_calltrace<T>::jacobian(const IR3& q) const {
  std::cout << "morphism::jacobian(q)\n";
  return T::jacobian(q);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_calltrace<T>::del_inverse(const IR3& q) const {
  std::cout << "morphism::del_inverse(q)\n";
  return T::del_inverse(q);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_calltrace<T>::tan_basis(const IR3& q) const {
  std::cout << "morphism::tan_basis(q)\n";
  return T::tan_basis(q);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_calltrace<T>::dual_basis(const IR3& q) const {
  std::cout << "morphism::dual_basis(q)\n";
  return T::dual_basis(q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_calltrace<T>::to_covariant(const IR3& A, const IR3& q) const {
  std::cout << "morphism::to_covariant(A, q)\n";
  return T::to_covariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_calltrace<T>::to_contravariant(const IR3& A, const IR3& q) const {
  std::cout << "morphism::to_contravariant(A, q)\n";
  return T::to_contravariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_calltrace<T>::from_covariant(const IR3& A, const IR3& q) const {
  std::cout << "morphism::from_covariant(A, q)\n";
  return T::from_covariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_calltrace<T>::from_contravariant(
    const IR3& A, const IR3& q) const {
  std::cout << "morphism::from_contravariant(A, q)\n";
  return T::from_contravariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_calltrace<T>::translation(const IR3& q, const IR3& delta) const {
  std::cout << "morphism::translation(q, delta)\n";
  return T::translation(q, delta);
}

}  // namespace gyronimo

#endif  // GYRONIMO_MORPHISM_CALLTRACE
