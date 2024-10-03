// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2024 Paulo Rodrigues.

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

#include <gyronimo/core/error.hh>
#include <gyronimo/core/generators.hh>

#include <gsl/gsl_multiroots.h>

#include <algorithm>
#include <functional>
#include <span>

namespace gyronimo {

//! Interface to GSL multiroot solvers without explicit derivatives.
/*!
    Provides a convenient interface for finding the zeros of continuously
    differentiable functions (or maps @f$f:\mathbb{R}^n \rightarrow
    \mathbb{R}^n@f$) via GSL routines. Allows arbitrary containers (constrained
    only to `SizedContiguousRange`) for storing data in @f$\mathbb{R}^n@f$ and
    as the arguments of the user-supplied functions whose root is to be found.
    Example of the intended usage:
    ```
    multiroot::settings_t my_settings = {
        .method = gsl_multiroot_fsolver_hybrids, .tolerance_abs = 1e-9,
        .tolerance_rel = 1e-9, .is_residual_tested = false, .iterations = 15};
    container_t guess = {1.0, ..., 2.2};
    container_t root =
        multiroot(my_settings)(my_root_f, guess);
    ```
    where `container_t` is any type following `SizedContiguousRange` and
    `my_root_f` is a `std::function<container_t(const container_t&)>` object
    containing the function whose root is being sought. Check the
    [GSL](https://www.gnu.org/software/gsl) documentation to better understand
    each method's properties and their eventual caveats.
*/
class multiroot {
 public:
  struct settings_t {
    const gsl_multiroot_fsolver_type* method;
    double tolerance_abs, tolerance_rel;
    bool is_residual_tested;
    size_t iterations;
  };
  multiroot(const settings_t& c) : settings_(c) {};
  const settings_t get_settings() const { return settings_; };

  template<SizedContiguousRange T_IRn>
  using user_map_t = typename std::function<T_IRn(const T_IRn&)>;
  template<SizedContiguousRange T_IRn>
  T_IRn operator()(user_map_t<T_IRn>& f, const T_IRn& guess) const;
 private:
  const settings_t settings_;

  using gsl_solver_t = gsl_multiroot_fsolver;
  inline bool is_converged(const gsl_solver_t* s) const;

  template<SizedContiguousRange T_IRn>
  static int translate_map_to_gsl(
      const gsl_vector* args_gsl, void* f_pointer, gsl_vector* eval_gsl);
  template<SizedContiguousRange T_IRn>
  auto allocate_gsl_objects(user_map_t<T_IRn>& f, const T_IRn& guess) const;
  inline void deallocate_gsl_objects(
      gsl_solver_t*, gsl_vector*, gsl_multiroot_function*) const;

  struct bad_alloc : public std::bad_alloc {};
  struct runtime_error : public std::runtime_error {
    runtime_error(const char* message) : std::runtime_error(message) {};
  };
};

template<SizedContiguousRange T_IRn>
T_IRn multiroot::operator()(user_map_t<T_IRn>& f, const T_IRn& guess) const {
  auto [solver, guess_gsl, struct_f_gsl] = allocate_gsl_objects(f, guess);
  for (auto iteration : std::views::iota(1u, settings_.iterations)) {
    int flag = gsl_multiroot_fsolver_iterate(solver);
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

inline bool multiroot::is_converged(const gsl_solver_t* solver) const {
  auto convergence_flag =
      (settings_.is_residual_tested ?
           gsl_multiroot_test_residual(solver->f, settings_.tolerance_abs) :
           gsl_multiroot_test_delta(
               solver->dx, solver->x, settings_.tolerance_abs,
               settings_.tolerance_rel));
  return (convergence_flag == GSL_SUCCESS ? true : false);
};

template<SizedContiguousRange T_IRn>
int multiroot::translate_map_to_gsl(
    const gsl_vector* args_gsl, void* f_pointer, gsl_vector* eval_gsl) {
  T_IRn args = generate_sized<T_IRn>(args_gsl->size);
  std::ranges::copy(std::span(args_gsl->data, args_gsl->size), args.begin());
  T_IRn eval = (*static_cast<user_map_t<T_IRn>*>(f_pointer))(args);
  std::ranges::copy(eval, eval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange T_IRn>
auto multiroot::allocate_gsl_objects(
    user_map_t<T_IRn>& f, const T_IRn& guess) const {
  const size_t n = guess.size();
  gsl_solver_t* solver = gsl_multiroot_fsolver_alloc(settings_.method, n);
  gsl_vector* guess_gsl = gsl_vector_alloc(n);
  auto* struct_f_gsl =
      new gsl_multiroot_function {&translate_map_to_gsl<T_IRn>, n, &f};
  if (!solver || !guess_gsl || !struct_f_gsl) throw bad_alloc();
  std::ranges::copy(guess, guess_gsl->data);
  gsl_multiroot_fsolver_set(solver, struct_f_gsl, guess_gsl);
  return std::tuple<gsl_solver_t*, gsl_vector*, gsl_multiroot_function*> {
      solver, guess_gsl, struct_f_gsl};
}

inline void multiroot::deallocate_gsl_objects(
    gsl_solver_t* solver, gsl_vector* guess_gsl,
    gsl_multiroot_function* struct_f_gsl) const {
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(guess_gsl);
  delete struct_f_gsl;
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_MULTIROOT
