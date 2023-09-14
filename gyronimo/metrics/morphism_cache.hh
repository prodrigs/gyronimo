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

// @morphism_cache.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_CACHE
#define GYRONIMO_MORPHISM_CACHE

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Cache layer over objects derived from `morphism`.
/*!
    Defines a depth-1 cache layer (or a decorator) for all member functions
    declared in the `morphism` interface. Because it derives from its template
    argument, it can be used to replace pointers or references to objects of
    that type in client code. The cache object is constructed with the same
    constructor call as its template argument.
*/
template<typename T> requires std::derived_from<T, morphism>
class morphism_cache : public T {
 public:
  template<typename... Args> morphism_cache(Args... args) : T(args...) {};
  virtual ~morphism_cache() {};
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
IR3 morphism_cache<T>::operator()(const IR3& q) const {
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::operator()(q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_cache<T>::inverse(const IR3& x) const {
  thread_local IR3 cache_x = un_init_, cache_eval = {0, 0, 0};
  if (x == cache_x) return cache_eval;
  cache_x = x;
  return cache_eval = T::inverse(x);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_cache<T>::del(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local dIR3 cache_eval = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::del(q);
}
template<typename T> requires std::derived_from<T, morphism>
ddIR3 morphism_cache<T>::ddel(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local ddIR3 cache_eval = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::ddel(q);
}
template<typename T> requires std::derived_from<T, morphism>
double morphism_cache<T>::jacobian(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local double cache_eval = 0;
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::jacobian(q);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_cache<T>::del_inverse(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local dIR3 cache_eval = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::del_inverse(q);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_cache<T>::tan_basis(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local dIR3 cache_eval = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::tan_basis(q);
}
template<typename T> requires std::derived_from<T, morphism>
dIR3 morphism_cache<T>::dual_basis(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local dIR3 cache_eval = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::dual_basis(q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_cache<T>::to_covariant(const IR3& A, const IR3& q) const {
  thread_local IR3 cache_q = un_init_, cache_A = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && A == cache_A) return cache_eval;
  cache_q = q;
  cache_A = A;
  return cache_eval = T::to_covariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_cache<T>::to_contravariant(const IR3& A, const IR3& q) const {
  thread_local IR3 cache_q = un_init_, cache_A = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && A == cache_A) return cache_eval;
  cache_q = q;
  cache_A = A;
  return cache_eval = T::to_contravariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_cache<T>::from_covariant(const IR3& A, const IR3& q) const {
  thread_local IR3 cache_q = un_init_, cache_A = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && A == cache_A) return cache_eval;
  cache_q = q;
  cache_A = A;
  return cache_eval = T::from_covariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_cache<T>::from_contravariant(const IR3& A, const IR3& q) const {
  thread_local IR3 cache_q = un_init_, cache_A = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && A == cache_A) return cache_eval;
  cache_q = q;
  cache_A = A;
  return cache_eval = T::from_contravariant(A, q);
}
template<typename T> requires std::derived_from<T, morphism>
IR3 morphism_cache<T>::translation(const IR3& q, const IR3& delta) const {
  thread_local IR3 cache_q = un_init_, cache_delta = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && delta == cache_delta) return cache_eval;
  cache_q = q;
  cache_delta = delta;
  return cache_eval = T::translation(q, delta);
}

}  // namespace gyronimo

#endif  // GYRONIMO_MORPHISM_CACHE
