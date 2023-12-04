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

// @metric_cache.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_CACHE
#define GYRONIMO_METRIC_CACHE

#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! Cache layer over objects derived from `metric_covariant`.
/*!
    Defines a depth-1 cache layer (or a decorator) for all member functions
    declared in the `metric_covariant` interface. Because it derives from its
    template argument, it can be used to replace pointers or references to
    objects of that type in client code. The cache object is constructed with
    the same constructor call as its template argument.
*/
template<typename T> requires std::derived_from<T, metric_covariant>
class metric_cache : public T {
 public:
  template<typename... Args> metric_cache(Args... args) : T(args...) {};
  virtual ~metric_cache() {};
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
 private:
  static inline IR3 un_init_ = {123456789., 987654321., 192837465.};
};

template<typename T> requires std::derived_from<T, metric_covariant>
SM3 metric_cache<T>::operator()(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local SM3 cache_eval = {0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::operator()(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
dSM3 metric_cache<T>::del(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local dSM3 cache_eval = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::del(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
double metric_cache<T>::jacobian(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local double cache_eval = 0;
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::jacobian(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_cache<T>::del_jacobian(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::del_jacobian(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
SM3 metric_cache<T>::inverse(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local SM3 cache_eval = {0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::inverse(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
dSM3 metric_cache<T>::del_inverse(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local dSM3 cache_eval = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::del_inverse(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
ddIR3 metric_cache<T>::christoffel_first_kind(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local ddIR3 cache_eval = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::christoffel_first_kind(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
ddIR3 metric_cache<T>::christoffel_second_kind(const IR3& q) const {
  thread_local IR3 cache_q = un_init_;
  thread_local ddIR3 cache_eval = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q) return cache_eval;
  cache_q = q;
  return cache_eval = T::christoffel_second_kind(q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_cache<T>::inertial_force(const IR3& q, const IR3& dot_q) const {
  thread_local IR3 cache_q = un_init_, cache_dot_q = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && dot_q == cache_dot_q) return cache_eval;
  cache_q = q;
  cache_dot_q = dot_q;
  return cache_eval = T::inertial_force(q, dot_q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_cache<T>::to_covariant(const IR3& B, const IR3& q) const {
  thread_local IR3 cache_q = un_init_, cache_B = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && B == cache_B) return cache_eval;
  cache_q = q;
  cache_B = B;
  return cache_eval = T::to_covariant(B, q);
}
template<typename T> requires std::derived_from<T, metric_covariant>
IR3 metric_cache<T>::to_contravariant(const IR3& B, const IR3& q) const {
  thread_local IR3 cache_q = un_init_, cache_B = un_init_;
  thread_local IR3 cache_eval = {0, 0, 0};
  if (q == cache_q && B == cache_B) return cache_eval;
  cache_q = q;
  cache_B = B;
  return cache_eval = T::to_contravariant(B, q);
}

}  // namespace gyronimo

#endif  // GYRONIMO_METRIC_CACHE
