// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues.

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

// @multiroot.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MULTIROOT
#define GYRONIMO_MULTIROOT

#include <span>
#include <ranges>
#include <algorithm>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gyronimo/core/error.hh>

namespace gyronimo {

//! Interface to [GSL](https://www.gnu.org/software/gsl) multiroot solver.
template<double Tol, size_t MaxIts>
class multiroot {
 public:
  template<typename UserArgs>
  using user_function_t = typename std::function<UserArgs(const UserArgs&)>;
  template<typename UserArgs>
  UserArgs operator()(user_function_t<UserArgs> f, const UserArgs& guess) const;
 private:
  template<typename UserArgs>
  static int translation_function(
      const gsl_vector *args_gsl, void *f, gsl_vector *eval_gsl);
  template<typename UserArgs>
  static std::tuple<gsl_multiroot_fsolver*, gsl_vector*> allocate_gsl_objects(
      user_function_t<UserArgs> f, const UserArgs& guess);
};

template<double Tol, size_t MaxIts>
template<typename UserArgs>
UserArgs multiroot<Tol, MaxIts>::operator()(
    std::function<UserArgs(const UserArgs&)> f, const UserArgs& guess) const {
  auto [solver, guess_gsl] = allocate_gsl_objects(f, guess);
  for(auto iteration : std::views::iota(1u, MaxIts)) {
    int flag = gsl_multiroot_fsolver_iterate(solver);
    switch(flag) {
      case GSL_ENOPROG: gyronimo::error(__func__, __FILE__, __LINE__,
          "iteration is stuck.", 1);
      case GSL_ENOPROGJ: gyronimo::error(__func__, __FILE__, __LINE__,
          "jacobian not improving the solution.", 1);
      case GSL_EBADFUNC: gyronimo::error(__func__, __FILE__, __LINE__,
          "singular user-supplied function (Inf/NaN).", 1);
    }
    if(gsl_multiroot_test_residual(solver->f, Tol) == GSL_SUCCESS) break;
  }
  if(gsl_multiroot_test_residual(solver->f, Tol) == GSL_CONTINUE)
      gyronimo::error(__func__, __FILE__, __LINE__,
          "still above tolerance after max iterations.", 1);
  UserArgs root = std::span<double>(solver->x->data, solver->x->size);
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(guess_gsl);
  return root;
}

template<double Tol, size_t MaxIts>
template<typename UserArgs>
int multiroot<Tol, MaxIts>::translation_function(
    const gsl_vector *args_gsl, void *f, gsl_vector *eval_gsl) {
  typedef std::function<UserArgs(const UserArgs&)> UserArgs_mapping;
  UserArgs args = std::span<double>(args_gsl->data, args_gsl->size);
  UserArgs eval = (*((user_function_t<UserArgs>*) f))(args);
  std::ranges::copy(eval, eval_gsl->data);
  return GSL_SUCCESS;
}

template<double Tol, size_t MaxIts>
template<typename UserArgs>
std::tuple<gsl_multiroot_fsolver*, gsl_vector*>
multiroot<Tol, MaxIts>::allocate_gsl_objects(
    user_function_t<UserArgs> f, const UserArgs& guess) {
  const size_t n = guess.size();
  const gsl_multiroot_fsolver_type *type = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(type, n);
  gsl_multiroot_function f_gsl = {&translation_function<UserArgs>, n, &f};
  gsl_vector *guess_gsl = gsl_vector_alloc(n);
  std::ranges::copy(guess, guess_gsl->data);
  gsl_multiroot_fsolver_set(solver, &f_gsl, guess_gsl);
  return {solver, guess_gsl};
}

} // end namespace gyronimo.

#endif // GYRONIMO_MULTIROOT
