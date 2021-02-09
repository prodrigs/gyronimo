// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @linear_combo.hh

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
  virtual ~linear_combo() {};
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
