// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2023 Manuel Assunção.

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

// @multiroot_c1.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MULTIROOT_C1
#define GYRONIMO_MULTIROOT_C1

#include <gyronimo/core/generators.hh>

#include <gsl/gsl_multiroots.h>

#include <algorithm>
#include <functional>
#include <span>

namespace gyronimo {

//! Interface to [GSL](https://www.gnu.org/software/gsl) C1 multiroot solver.
/*!
    Extends `multiroot` functionality to continuously differentiable functions
    with user-supplied derivatives. Examples of intended usage:

    - Most general case, may not be optimal when evaluating derivatives:
    ```
    container_t guess = {1.0, ..., 2.2};
    container_t root = multiroot_c1(gsl_multiroot_fdfsolver_hybridj, 1e-9, 75)(
        my_root_f, my_root_df, my_root_fdf, guess);
    ```
    - Simplified case, recommended only when the computation of the function
    cannot be recycled for the computation of the jacobian:
    ```
    container_t guess = {1.0, ..., 2.2};
    container_t root = multiroot_c1(gsl_multiroot_fdfsolver_gnewton, 1e-9, 75)(
        my_root_f, my_root_df, guess);
    ```
    - Simplest case, works only with the `newton` method:
    ```
    container_t guess = {1.0, ..., 2.2};
    container_t root = multiroot_c1(gsl_multiroot_fdfsolver_newton, 1e-9, 75)(
        my_root_fdf, guess);
    ```
    Here, `container_t` and `container_d_t` are any storage types following
    `SizedContiguousRange`, `my_root_f` is a `std::function<container_t(const
    container_t&)>` object holding the function to find the root of,
    `my_root_df` is a `std::function<container_d_t(const container_t&)>` object
    holding the jacobian of the function to find the root of, and `my_root_fdf`
    is a `std::function<std::pair<container_t,container_d_t>(const
    container_t&)>` object containing both the function to find the root of and
    its jacobian. Each `multiroot_c1::solver` attribute is based on a GSL
    fdfsolver_type. Check the GSL documentation to better understand each
    method's properties, and their eventual caveats.
*/
class multiroot_c1 {
 public:
  struct runtime_error : public std::runtime_error {
    runtime_error(const char* message) : std::runtime_error(message) {};
  };
  struct configuration_t {
    const gsl_multiroot_fdfsolver_type* method;
    const double tolerance_abs, tolerance_rel;
    const bool is_residual_tested;
    const size_t iterations;
  };

  multiroot_c1(const configuration_t& c) : configuration_(c) {};

  template<SizedContiguousRange UserArgs>
  using user_function_t = typename std::function<UserArgs(const UserArgs&)>;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  using user_dfunction_t = typename std::function<UserDArgs(const UserArgs&)>;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  using user_fdfunction_t =
      typename std::function<std::pair<UserArgs, UserDArgs>(const UserArgs&)>;

  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  UserArgs operator()(
      user_fdfunction_t<UserArgs, UserDArgs>& fdf, const UserArgs& guess) const;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  UserArgs operator()(
      user_function_t<UserArgs>& f, user_dfunction_t<UserArgs, UserDArgs>& df,
      const UserArgs& guess) const;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  UserArgs operator()(
      user_function_t<UserArgs>& f, user_dfunction_t<UserArgs, UserDArgs>& df,
      user_fdfunction_t<UserArgs, UserDArgs>& fdf, const UserArgs& guess) const;

  const configuration_t get_configuration() const { return configuration_; };
 private:
  const configuration_t configuration_;

  inline bool is_iteration_converged(const gsl_multiroot_fdfsolver* s) const;
  inline void deallocate_gsl_objects(
      gsl_multiroot_fdfsolver* solver, gsl_vector* guess_gsl,
      gsl_multiroot_function_fdf* struct_fdf_gsl) const;

  template<SizedContiguousRange UserArgs>
  UserArgs solve_kernel(gsl_multiroot_fdfsolver* solver) const;

  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  struct function_pack_t {
    user_function_t<UserArgs>* f;
    user_dfunction_t<UserArgs, UserDArgs>* df;
    user_fdfunction_t<UserArgs, UserDArgs>* fdf;
  };
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  auto allocate_gsl_objects(
      function_pack_t<UserArgs, UserDArgs>& fpack, const UserArgs& guess) const;

  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  static int translation_function(
      const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl);
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  static int translation_dfunction(
      const gsl_vector* args_gsl, void* fpack, gsl_matrix* deval_gsl);
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  static int translation_fdfunction(
      const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl,
      gsl_matrix* deval_gsl);
};

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
UserArgs multiroot_c1::operator()(
    user_function_t<UserArgs>& f, user_dfunction_t<UserArgs, UserDArgs>& df,
    user_fdfunction_t<UserArgs, UserDArgs>& fdf, const UserArgs& guess) const {
  function_pack_t<UserArgs, UserDArgs> pack = {&f, &df, &fdf};
  auto [solver, guess_gsl, struct_fdf_gsl] =
      allocate_gsl_objects<UserArgs, UserDArgs>(pack, guess);
  UserArgs root = solve_kernel<UserArgs>(solver);
  deallocate_gsl_objects(solver, guess_gsl, struct_fdf_gsl);
  return root;
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
UserArgs multiroot_c1::operator()(
    user_function_t<UserArgs>& f, user_dfunction_t<UserArgs, UserDArgs>& df,
    const UserArgs& guess) const {
  user_fdfunction_t<UserArgs, UserDArgs> fdf =
      [&](UserArgs args) -> std::pair<UserArgs, UserDArgs> {
    return {f(args), df(args)};
  };
  return (*this)(f, df, fdf, guess);
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
UserArgs multiroot_c1::operator()(
    user_fdfunction_t<UserArgs, UserDArgs>& fdf, const UserArgs& guess) const {
  user_function_t<UserArgs> f = [&](UserArgs args) -> UserArgs {
    std::pair<UserArgs, UserDArgs> eval = fdf(args);
    return eval.first;
  };
  user_dfunction_t<UserArgs, UserDArgs> df = [&](UserArgs args) -> UserDArgs {
    std::pair<UserArgs, UserDArgs> eval = fdf(args);
    return eval.second;
  };
  return (*this)(f, df, fdf, guess);
}

template<SizedContiguousRange UserArgs>
UserArgs multiroot_c1::solve_kernel(gsl_multiroot_fdfsolver* solver) const {
  for (auto iteration : std::views::iota(1u, configuration_.iterations)) {
    int flag = gsl_multiroot_fdfsolver_iterate(solver);
    switch (flag) {
      case GSL_ENOPROG: throw runtime_error("iteration is stuck.");
      case GSL_ENOPROGJ: throw runtime_error("jacobian not improving.");
      case GSL_EBADFUNC: throw runtime_error("singular function (Inf/NaN).");
    }
    if (this->is_iteration_converged(solver)) break;
  }
  if (!this->is_iteration_converged(solver))
    throw runtime_error("max iterations exceeded.");

  UserArgs root = generate_sized<UserArgs>(solver->x->size);
  std::ranges::copy(std::span(solver->x->data, solver->x->size), root.begin());
  return root;
}

inline bool multiroot_c1::is_iteration_converged(const gsl_multiroot_fdfsolver* solver) const {
  auto convergence_flag = (configuration_.is_residual_tested ?
      gsl_multiroot_test_residual(solver->f, configuration_.tolerance_abs) :
      gsl_multiroot_test_delta(solver->dx, solver->x, configuration_.tolerance_abs, configuration_.tolerance_rel));
  return (convergence_flag == GSL_SUCCESS ? true : false);
};

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
int multiroot_c1::translation_function(
    const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl) {
  UserArgs args = generate_sized<UserArgs>(args_gsl->size);
  std::ranges::copy(std::span(args_gsl->data, args_gsl->size), args.begin());
  auto* pk = static_cast<function_pack_t<UserArgs, UserDArgs>*>(fpack);
  UserArgs eval = (*pk->f)(args);
  std::ranges::copy(eval, eval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
int multiroot_c1::translation_dfunction(
    const gsl_vector* args_gsl, void* fpack, gsl_matrix* deval_gsl) {
  UserArgs args = generate_sized<UserArgs>(args_gsl->size);
  std::ranges::copy(std::span(args_gsl->data, args_gsl->size), args.begin());
  auto* pk = static_cast<function_pack_t<UserArgs, UserDArgs>*>(fpack);
  UserDArgs deval = (*pk->df)(args);
  std::ranges::copy(deval, deval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
int multiroot_c1::translation_fdfunction(
    const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl,
    gsl_matrix* deval_gsl) {
  UserArgs args = generate_sized<UserArgs>(args_gsl->size);
  std::ranges::copy(std::span(args_gsl->data, args_gsl->size), args.begin());
  auto* pk = static_cast<function_pack_t<UserArgs, UserDArgs>*>(fpack);
  std::pair<UserArgs, UserDArgs> eval = (*pk->fdf)(args);
  std::ranges::copy(eval.first, eval_gsl->data);
  std::ranges::copy(eval.second, deval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
auto multiroot_c1::allocate_gsl_objects(
    function_pack_t<UserArgs, UserDArgs>& fpack, const UserArgs& guess) const {
  const size_t n = guess.size();
  gsl_multiroot_fdfsolver* solver = gsl_multiroot_fdfsolver_alloc(configuration_.method, n);
  gsl_vector* guess_gsl = gsl_vector_alloc(n);
  auto* struct_fdf_gsl = new gsl_multiroot_function_fdf {
      &translation_function<UserArgs, UserDArgs>,
      &translation_dfunction<UserArgs, UserDArgs>,
      &translation_fdfunction<UserArgs, UserDArgs>, n, &fpack};
  if (!solver || !guess_gsl || !struct_fdf_gsl)
    error(__func__, __FILE__, __LINE__, "gsl object allocation failed.", 1);
  std::ranges::copy(guess, guess_gsl->data);
  gsl_multiroot_fdfsolver_set(solver, struct_fdf_gsl, guess_gsl);
  return std::tuple<
      gsl_multiroot_fdfsolver*, gsl_vector*, gsl_multiroot_function_fdf*> {
      solver, guess_gsl, struct_fdf_gsl};
}

inline void multiroot_c1::deallocate_gsl_objects(
    gsl_multiroot_fdfsolver* solver, gsl_vector* guess_gsl,
    gsl_multiroot_function_fdf* struct_fdf_gsl) const {
  gsl_multiroot_fdfsolver_free(solver);
  gsl_vector_free(guess_gsl);
  delete struct_fdf_gsl;
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_MULTIROOT_C1
