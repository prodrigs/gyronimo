// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Manuel Assunção and Paulo Rodrigues.

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

// @metric_connected.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_CONNECTED
#define GYRONIMO_METRIC_CONNECTED

#include <gyronimo/core/error.hh>
#include <gyronimo/metrics/metric_covariant.hh>
#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! Covariant metric connected to a defining `gyronimo::morphism`.
/*!
    Unlike `metric_covariant`, here the underlying transformation from the
    specific curvilinear coordinates into the cartesian space is established by
    the `gyronimo::morphism` object supplied at the construction. This
    information is sufficient to override every abstract member function of the
    parent class (and some non-abstract member functions also) with
    general-purpose implementations. In this sense, `metric_connected` is
    already a full functional (i.e., instantiable) object. However, it can be
    used as a parent for derived classes intended to override some of its member
    functions in order to have them specialised (and thus optimised) according
    to the particular properties enjoyed by specific coordinate sets.
*/
class metric_connected : public metric_covariant {
 public:
  metric_connected(const morphism* m) : my_morphism_(m) {};
  virtual ~metric_connected() override {};
  virtual SM3 operator()(const IR3& q) const override;
  virtual dSM3 del(const IR3& q) const override;
  virtual double jacobian(const IR3& q) const override;
  virtual IR3 del_jacobian(const IR3& q) const override;
  virtual ddIR3 christoffel_first_kind(const IR3& q) const override;
  virtual ddIR3 christoffel_second_kind(const IR3& q) const override;

  const morphism* my_morphism() const { return my_morphism_; };
 private:
  const morphism* my_morphism_;
};

//! General-purpose jacobian, as inherited from parent `morphism`.
inline double metric_connected::jacobian(const IR3& q) const {
  return my_morphism_->jacobian(q);
}

//! Christoffel @f$\Gamma_{kij}=\mathbf{e}_k\cdot\partial^2_{ij}\mathbf{x}@f$.
inline ddIR3 metric_connected::christoffel_first_kind(const IR3& q) const {
  return contraction<first>(my_morphism_->del(q), my_morphism_->ddel(q));
}

//! Christoffel @f$\Gamma^k_{ij}=\mathbf{e}^k\cdot\partial^2_{ij}\mathbf{x}@f$.
inline ddIR3 metric_connected::christoffel_second_kind(const IR3& q) const {
  return contraction<second>(
      my_morphism_->del_inverse(q), my_morphism_->ddel(q));
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_METRIC_CONNECTED
