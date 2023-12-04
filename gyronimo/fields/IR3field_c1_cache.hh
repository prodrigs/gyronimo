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

// @IR3field_c1_cache.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_IR3FIELD_C1_CACHE
#define GYRONIMO_IR3FIELD_C1_CACHE

#include <gyronimo/fields/IR3field_c1.hh>

namespace gyronimo {

//! Calltrace layer over objects derived from `IR3field_c1`.
/*!
    Defines a depth-1 cache layer (or a decorator) for all member functions
    declared in the `IR3field_c1` interface. Because it derives from its
    template argument, it can be used to replace pointers or references to
    objects of that type in client code. The cache object is constructed with
    the same constructor call as its template argument.
*/
template<typename T> requires std::derived_from<T, IR3field_c1>
class IR3field_c1_cache : public T {
 public:
  template<typename... Args> IR3field_c1_cache(Args... args) : T(args...) {};
  virtual ~IR3field_c1_cache() {};
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
 private:
  static inline IR3 un_init_ = {123456789., 987654321., 192837465.};
};

template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::contravariant(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::contravariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::covariant(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::covariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
double IR3field_c1_cache<T>::magnitude(const IR3& q, double t) const {
  thread_local double cache_t = -1, cache_eval = 0;
  thread_local IR3 cache_q = un_init_;
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::magnitude(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::covariant_versor(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::covariant_versor(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::contravariant_versor(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::contravariant_versor(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
dIR3 IR3field_c1_cache<T>::del_contravariant(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_;
  thread_local dIR3 cache_eval = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::del_contravariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::partial_t_contravariant(
    const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::partial_t_contravariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::del_magnitude(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::del_magnitude(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
double IR3field_c1_cache<T>::partial_t_magnitude(const IR3& q, double t) const {
  thread_local double cache_t = -1, cache_eval = 0;
  thread_local IR3 cache_q = un_init_;
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::partial_t_magnitude(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
dIR3 IR3field_c1_cache<T>::del_covariant(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_;
  thread_local dIR3 cache_eval = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::del_covariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::partial_t_covariant(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::partial_t_covariant(q, t);
}
template<typename T> requires std::derived_from<T, IR3field_c1>
IR3 IR3field_c1_cache<T>::curl(const IR3& q, double t) const {
  thread_local double cache_t = -1;
  thread_local IR3 cache_q = un_init_, cache_eval = {0, 0, 0};
  if (q == cache_q && t == cache_t) return cache_eval;
  cache_q = q;
  cache_t = t;
  return cache_eval = T::curl(q, t);
}

}  // namespace gyronimo

#endif  // GYRONIMO_IR3FIELD_C1_CACHE
