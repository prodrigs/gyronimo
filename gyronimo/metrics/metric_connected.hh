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

//! Covariant metric connected to a defining `morphism`.
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

//! General Christoffel symbol @f$\Gamma_{kij}@f$ from parent `morphism`.
/*!
    @f{equation*}{\Gamma_{kij} = \textbf{e}_k \cdot
          \frac{\partial^2 \textbf{x}}{\partial q^i \, \partial q^j}@f}
*/
inline ddIR3 metric_connected::christoffel_first_kind(const IR3& q) const {
  return contraction<first, first>(my_morphism_->del(q), my_morphism_->ddel(q));
}

//! General Christoffel symbol @f$\Gamma^k_{ij}@f$ from parent `morphism`.
/*!
    @f{equation*}{\Gamma^k_{ij} = \textbf{e}^k \cdot
          \frac{\partial^2 \textbf{x}}{\partial q^i \, \partial q^j}@f}
*/
inline ddIR3 metric_connected::christoffel_second_kind(const IR3& q) const {
  return contraction<second, first>(
      my_morphism_->del_inverse(q), my_morphism_->ddel(q));
}

}  // end namespace gyronimo.

#endif  // GYRONIMO_METRIC_CONNECTED
