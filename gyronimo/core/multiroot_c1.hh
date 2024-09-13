// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2023-2024 Manuel Assunção and Paulo Rodrigues.

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

//! Interface to GSL multiroot solvers with explicit derivatives supplied.
/*!
    Extends `multiroot` functionality to continuously differentiable functions
    (or maps @f$f:\mathbb{R}^n \rightarrow \mathbb{R}^n@f$) with user-supplied
    derivatives. Some calling sequences allow the simultaneous evaluations of
    functions and derivatives for efficiency. Examples of intended usage:

    1. Simplest case, may be inefficient for methods other than `newton`:
    ```
    multiroot_c1::settings_t my_settings = {
        .method = gsl_multiroot_fdfsolver_newton, .tolerance_abs = 1e-9,
        .tolerance_rel = 1e-9, .is_residual_tested = false, .iterations = 15};
    container_t my_guess = {1.0, ..., 2.2};
    container_t root = multiroot_c1(my_settings)(my_root_fdf, my_guess);
    ```
    2. General case, may be more efficient for methods other than `newton`,
    brings no additional benefits for the latter:
    ```
    container_t root = multiroot_c1(my_settings)(
        my_root_f, my_root_df, my_root_fdf, my_guess);
    ```
    3. Simplified case, crafts a function-derivative combination from supplied
    my_root_f and my_root_df for the sake of lazy users:
    ```
    container_t root =
        multiroot_c1(my_settings)(my_root_f, my_root_df, my_guess);
    ```

    Here, `container_t` and `container_d_t` are any storage types following
    `SizedContiguousRange`, `my_root_f` is a `std::function<container_t(const
    container_t&)>` object holding the map to find the root of, `my_root_df` is
    a `std::function<container_d_t(const container_t&)>` object holding the
    jacobian of the map to find the root of, and `my_root_fdf` is a
    `std::function<std::pair<container_t,container_d_t>(const container_t&)>`
    object containing both the map and its jacobian (computed together for
    efficiency, if any). The attribute `multiroot_c1::settings_t::method` can
    hold any valid GSL method (i.e., any defined pointer
    `gsl_multiroot_fdfsolver_type*`). Check the
    [GSL](https://www.gnu.org/software/gsl) documentation to better understand
    each method's properties and their eventual caveats.
*/
class multiroot_c1 {
 public:
  struct settings_t {
    const gsl_multiroot_fdfsolver_type* method;
    const double tolerance_abs, tolerance_rel;
    const bool is_residual_tested;
    const size_t iterations;
  };
  multiroot_c1(const settings_t& c) : settings_(c) {};
  const settings_t get_settings() const { return settings_; };

  template<SizedContiguousRange T_IRn>
  using user_map_t = typename std::function<T_IRn(const T_IRn&)>;
  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  using user_dmap_t = typename std::function<T_dIRn(const T_IRn&)>;
  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  using user_combo_t =
      typename std::function<std::pair<T_IRn, T_dIRn>(const T_IRn&)>;

  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  T_IRn operator()(
      const user_combo_t<T_IRn, T_dIRn>& fdf, const T_IRn& guess) const;
  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  T_IRn operator()(
      const user_map_t<T_IRn>& f, const user_dmap_t<T_IRn, T_dIRn>& df,
      const T_IRn& guess) const;
  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  T_IRn operator()(
      const user_map_t<T_IRn>& f, const user_dmap_t<T_IRn, T_dIRn>& df,
      const user_combo_t<T_IRn, T_dIRn>& fdf, const T_IRn& guess) const;
 private:
  const settings_t settings_;

  using gsl_solver_t = gsl_multiroot_fdfsolver;
  inline bool is_converged(const gsl_solver_t* s) const;
  inline void deallocate_gsl_objects(
      gsl_solver_t* solver, gsl_vector* guess_gsl,
      gsl_multiroot_function_fdf* struct_fdf_gsl) const;

  template<SizedContiguousRange T_IRn>
  T_IRn solve_kernel(gsl_solver_t* solver) const;

  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  struct map_pack_t {
    const user_map_t<T_IRn>* f;
    const user_dmap_t<T_IRn, T_dIRn>* df;
    const user_combo_t<T_IRn, T_dIRn>* fdf;
  };
  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  auto allocate_gsl_objects(
      map_pack_t<T_IRn, T_dIRn>& map_pack, const T_IRn& guess) const;

  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  static int translate_map_to_gsl(
      const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl);
  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  static int translate_dmap_to_gsl(
      const gsl_vector* args_gsl, void* fpack, gsl_matrix* deval_gsl);
  template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
  static int translate_combo_to_gsl(
      const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl,
      gsl_matrix* deval_gsl);

  struct bad_alloc : public std::bad_alloc {};
  struct runtime_error : public std::runtime_error {
    runtime_error(const char* message) : std::runtime_error(message) {};
  };
};

template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
T_IRn multiroot_c1::operator()(
    const user_map_t<T_IRn>& f, const user_dmap_t<T_IRn, T_dIRn>& df,
    const user_combo_t<T_IRn, T_dIRn>& fdf, const T_IRn& guess) const {
  map_pack_t<T_IRn, T_dIRn> map_pack = {&f, &df, &fdf};
  auto [solver, guess_in_gsl, struct_fdf_gsl] =
      allocate_gsl_objects<T_IRn, T_dIRn>(map_pack, guess);
  T_IRn root = this->solve_kernel<T_IRn>(solver);
  deallocate_gsl_objects(solver, guess_in_gsl, struct_fdf_gsl);
  return root;
}

template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
T_IRn multiroot_c1::operator()(
    const user_map_t<T_IRn>& f, const user_dmap_t<T_IRn, T_dIRn>& df,
    const T_IRn& guess) const {
  user_combo_t<T_IRn, T_dIRn> crafted_fdf = [&](const T_IRn& args) {
    return {f(args), df(args)};
  };
  return (*this)(f, df, crafted_fdf, guess);
}

template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
T_IRn multiroot_c1::operator()(
    const user_combo_t<T_IRn, T_dIRn>& fdf, const T_IRn& guess) const {
  user_map_t<T_IRn> crafted_f = [&](const T_IRn& args) {
    return fdf(args).first;
  };
  user_dmap_t<T_IRn, T_dIRn> crafted_df = [&](const T_IRn& args) {
    return fdf(args).second;
  };
  return (*this)(crafted_f, crafted_df, fdf, guess);
}

template<SizedContiguousRange T_IRn>
T_IRn multiroot_c1::solve_kernel(gsl_solver_t* solver) const {
  for (auto iteration : std::views::iota(1u, settings_.iterations)) {
    int flag = gsl_multiroot_fdfsolver_iterate(solver);
    switch (flag) {
      case GSL_ENOPROG: throw runtime_error("iteration is stuck.");
      case GSL_ENOPROGJ: throw runtime_error("jacobian not improving.");
      case GSL_EBADFUNC: throw runtime_error("singular function (Inf/NaN).");
    }
    if (this->is_converged(solver)) break;
  }
  if (!this->is_converged(solver)) throw runtime_error("iterations exceeded.");
  T_IRn root = generate_sized<T_IRn>(solver->x->size);
  std::ranges::copy(std::span(solver->x->data, solver->x->size), root.begin());
  return root;
}

inline bool multiroot_c1::is_converged(const gsl_solver_t* solver) const {
  auto convergence_flag =
      (settings_.is_residual_tested ?
           gsl_multiroot_test_residual(solver->f, settings_.tolerance_abs) :
           gsl_multiroot_test_delta(
               solver->dx, solver->x, settings_.tolerance_abs,
               settings_.tolerance_rel));
  return (convergence_flag == GSL_SUCCESS ? true : false);
};

template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
int multiroot_c1::translate_map_to_gsl(
    const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl) {
  T_IRn args = generate_sized<T_IRn>(args_gsl->size);
  std::ranges::copy(std::span(args_gsl->data, args_gsl->size), args.begin());
  auto* pk = static_cast<map_pack_t<T_IRn, T_dIRn>*>(fpack);
  T_IRn eval = (*pk->f)(args);
  std::ranges::copy(eval, eval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
int multiroot_c1::translate_dmap_to_gsl(
    const gsl_vector* args_gsl, void* fpack, gsl_matrix* deval_gsl) {
  T_IRn args = generate_sized<T_IRn>(args_gsl->size);
  std::ranges::copy(std::span(args_gsl->data, args_gsl->size), args.begin());
  auto* pk = static_cast<map_pack_t<T_IRn, T_dIRn>*>(fpack);
  T_dIRn deval = (*pk->df)(args);
  std::ranges::copy(deval, deval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
int multiroot_c1::translate_combo_to_gsl(
    const gsl_vector* args_gsl, void* fpack, gsl_vector* eval_gsl,
    gsl_matrix* deval_gsl) {
  T_IRn args = generate_sized<T_IRn>(args_gsl->size);
  std::ranges::copy(std::span(args_gsl->data, args_gsl->size), args.begin());
  auto* pk = static_cast<map_pack_t<T_IRn, T_dIRn>*>(fpack);
  std::pair<T_IRn, T_dIRn> eval = (*pk->fdf)(args);
  std::ranges::copy(eval.first, eval_gsl->data);
  std::ranges::copy(eval.second, deval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange T_IRn, SizedContiguousRange T_dIRn>
auto multiroot_c1::allocate_gsl_objects(
    map_pack_t<T_IRn, T_dIRn>& map_pack, const T_IRn& guess) const {
  const size_t n = guess.size();
  gsl_solver_t* solver = gsl_multiroot_fdfsolver_alloc(settings_.method, n);
  gsl_vector* guess_gsl = gsl_vector_alloc(n);
  auto* struct_fdf_gsl = new gsl_multiroot_function_fdf {
      &translate_map_to_gsl<T_IRn, T_dIRn>,
      &translate_dmap_to_gsl<T_IRn, T_dIRn>,
      &translate_combo_to_gsl<T_IRn, T_dIRn>, n, &map_pack};
  if (!solver || !guess_gsl || !struct_fdf_gsl) throw bad_alloc();
  std::ranges::copy(guess, guess_gsl->data);
  gsl_multiroot_fdfsolver_set(solver, struct_fdf_gsl, guess_gsl);
  return std::tuple<gsl_solver_t*, gsl_vector*, gsl_multiroot_function_fdf*> {
      solver, guess_gsl, struct_fdf_gsl};
}

inline void multiroot_c1::deallocate_gsl_objects(
    gsl_solver_t* solver, gsl_vector* guess_gsl,
    gsl_multiroot_function_fdf* struct_fdf_gsl) const {
  gsl_multiroot_fdfsolver_free(solver);
  gsl_vector_free(guess_gsl);
  delete struct_fdf_gsl;
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_MULTIROOT_C1
