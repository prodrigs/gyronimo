// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @linear_combo.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_LINEAR_COMBO
#define GYRONIMO_LINEAR_COMBO

#include <gyronimo/core/error.hh>
#include <gyronimo/fields/IR3field.hh>

namespace gyronimo {

//! Linear combination of `N` IR3 fields (*shared* coordinates).
template<size_t N>
class linear_combo : public IR3field {
 public:
  linear_combo(
      const std::array<const IR3field*, N>& p,
      const metric_covariant* g, double m_factor, double t_factor);
  virtual ~linear_combo() override {};
  virtual IR3 contravariant(const IR3& position, double time) const override;
 private:
  std::array<const IR3field*, N> field_set_;
  std::array<double, N> m_ratio_, t_ratio_;
};

template<size_t N>
linear_combo<N>::linear_combo(
    const std::array<const IR3field*, N>& p,
    const metric_covariant* g, double m_factor, double t_factor)
    : IR3field(m_factor, t_factor, g), field_set_(p) {
  for (size_t i = 0; i < N; i++) {
    if (field_set_[i]->metric() != g)
        error(__func__, __FILE__, __LINE__, "incompatible metrics.", 1);
    m_ratio_[i] = field_set_[i]->m_factor()/m_factor;
    t_ratio_[i] = t_factor/field_set_[i]->t_factor();
  }
}
template<size_t N>
IR3 linear_combo<N>::contravariant(const IR3& position, double time) const {
  IR3 acc = {0, 0, 0};
  for (size_t i = 0; i < N; i++) acc +=
      m_ratio_[i]*field_set_[i]->contravariant(position, t_ratio_[i]*time);
  return acc;
}

} // end namespace gyronimo.

#endif // GYRONIMO_LINEAR_COMBO
