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

#include <ranges>
#include <algorithm>
#include <gsl/gsl_multiroots.h>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/generators.hh>

namespace gyronimo {

//! Interface to [GSL](https://www.gnu.org/software/gsl) multiroot solver.
/*!
    Intended usage:
    ```
    container_t guess = {1.0, ..., 2.2};
    container_t root = multiroot(1.0e-09, 75)(my_rooting_function, guess);
    ```
    where `container_t` is any type following `SizedContiguousRange` and
    `my_rooting_function` is a `std::function<container_t(container_t)>` object
    containing the function to find the root of.
*/
class multiroot {
 public:
  template<SizedContiguousRange UserArgs>
  using user_function_t = typename std::function<UserArgs(const UserArgs&)>;
  multiroot(double tolerance, size_t max_iterations)
      : max_iterations_(max_iterations), tolerance_(tolerance) {};
  template<SizedContiguousRange UserArgs>
  UserArgs operator()(
      user_function_t<UserArgs>& f, const UserArgs& guess) const;
 private:
  size_t max_iterations_;
  double tolerance_;
  template<SizedContiguousRange UserArgs>
  static int translation_function(
      const gsl_vector *args_gsl, void *f_pointer, gsl_vector *eval_gsl);
  template<SizedContiguousRange UserArgs>
  auto allocate_gsl_objects(
      user_function_t<UserArgs>& f, const UserArgs& guess) const;
  void deallocate_gsl_objects(
      gsl_multiroot_fsolver*, gsl_vector*, gsl_multiroot_function*) const;
};

template<SizedContiguousRange UserArgs>
UserArgs multiroot::operator()(
    user_function_t<UserArgs>& f, const UserArgs& guess) const {
  auto [solver, guess_gsl, struct_f_gsl] = allocate_gsl_objects(f, guess);
  for(auto iteration : std::views::iota(1u, max_iterations_)) {
    int flag = gsl_multiroot_fsolver_iterate(solver);
    switch(flag) {
      case GSL_ENOPROG: error(__func__, __FILE__, __LINE__,
          "iteration is stuck.", 1);
      case GSL_ENOPROGJ: error(__func__, __FILE__, __LINE__,
          "jacobian not improving the solution.", 1);
      case GSL_EBADFUNC: error(__func__, __FILE__, __LINE__,
          "singular user-supplied function (Inf/NaN).", 1);
    }
    if(gsl_multiroot_test_residual(solver->f, tolerance_) == GSL_SUCCESS) break;
  }
  if(gsl_multiroot_test_residual(solver->f, tolerance_) == GSL_CONTINUE)
      error(__func__, __FILE__, __LINE__,
          "still above tolerance after max iterations.", 1);
  UserArgs root = generate_sized<UserArgs>(solver->x->size);
  std::copy(solver->x->data, solver->x->data + solver->x->size, root.begin());
  deallocate_gsl_objects(solver, guess_gsl, struct_f_gsl);
  return root;
}

template<SizedContiguousRange UserArgs>
int multiroot::translation_function(
    const gsl_vector *args_gsl, void *f_pointer, gsl_vector *eval_gsl) {
  UserArgs args = generate_sized<UserArgs>(args_gsl->size);
  std::copy(args_gsl->data, args_gsl->data + args_gsl->size, args.begin());
  UserArgs eval = (*static_cast<user_function_t<UserArgs>*>(f_pointer))(args);
  std::ranges::copy(eval, eval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange UserArgs>
auto multiroot::allocate_gsl_objects(
    user_function_t<UserArgs>& f, const UserArgs& guess) const {
  const size_t n = guess.size();
  const gsl_multiroot_fsolver_type *type = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(type, n);
  gsl_vector *guess_gsl = gsl_vector_alloc(n);
  auto *struct_f_gsl =
      new gsl_multiroot_function {&translation_function<UserArgs>, n, &f};
  if(!solver || !guess_gsl || !struct_f_gsl) error(
      __func__, __FILE__, __LINE__, "gsl object allocation failed.", 1);
  std::ranges::copy(guess, guess_gsl->data);
  gsl_multiroot_fsolver_set(solver, struct_f_gsl, guess_gsl);
  return std::tuple<
      gsl_multiroot_fsolver*, gsl_vector*, gsl_multiroot_function*
          >({solver, guess_gsl, struct_f_gsl});
}

void multiroot::deallocate_gsl_objects(
    gsl_multiroot_fsolver* solver,
    gsl_vector* guess_gsl, gsl_multiroot_function* struct_f_gsl) const {
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(guess_gsl);
  delete struct_f_gsl;
}

} // end namespace gyronimo.

#endif // GYRONIMO_MULTIROOT
