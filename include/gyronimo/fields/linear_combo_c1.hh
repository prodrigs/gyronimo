// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @linear_combo_c1.hh

#ifndef GYRONIMO_LINEAR_COMBO_C1
#define GYRONIMO_LINEAR_COMBO_C1

#include <gyronimo/core/error.hh>
#include <gyronimo/fields/IR3field_c1.hh>

namespace gyronimo {

//! Linear combination of `N` differentiable IR3 fields (*shared* coordinates).
template<size_t N>
class linear_combo_c1 : public IR3field_c1 {
 public:
  linear_combo_c1(
      const std::array<const IR3field_c1*, N>& p,
      const metric_covariant* g, double m_factor, double t_factor);
  virtual ~linear_combo_c1() override {};
  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override;
 private:
  std::array<const IR3field_c1*, N> field_set_;
  std::array<double, N> m_ratio_, t_ratio_;
};

template<size_t N>
linear_combo_c1<N>::linear_combo_c1(
    const std::array<const IR3field_c1*, N>& p,
    const metric_covariant* g, double m_factor, double t_factor)
    : IR3field_c1(m_factor, t_factor, g), field_set_(p) {
  for (size_t i = 0; i < N; i++) {
    if (field_set_[i]->metric() != g)
        error(__func__, __FILE__, __LINE__, "incompatible metrics.", 1);
    m_ratio_[i] = field_set_[i]->m_factor()/m_factor;
    t_ratio_[i] = t_factor/field_set_[i]->t_factor();
  }
}
template<size_t N>
IR3 linear_combo_c1<N>::contravariant(const IR3& position, double time) const {
  IR3 acc = {0, 0, 0};
  for (size_t i = 0; i < N; i++) acc +=
      m_ratio_[i]*field_set_[i]->contravariant(position, t_ratio_[i]*time);
  return acc;
}
template<size_t N>
IR3 linear_combo_c1<N>::del_contravariant(
    const IR3& position, double time) const {
  dIR3 acc = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (size_t i = 0; i < N; i++) {
    dIR3 del_i = field_set_[i]->del_contravariant(position, t_ratio_[i]*time);
    acc[dIR3::uu] += m_ratio_[i]*del_i[dIR3::uu];
    acc[dIR3::uv] += m_ratio_[i]*del_i[dIR3::uv];
    acc[dIR3::uw] += m_ratio_[i]*del_i[dIR3::uw];
    acc[dIR3::vu] += m_ratio_[i]*del_i[dIR3::vu];
    acc[dIR3::vv] += m_ratio_[i]*del_i[dIR3::vv];
    acc[dIR3::vw] += m_ratio_[i]*del_i[dIR3::vw];
    acc[dIR3::wu] += m_ratio_[i]*del_i[dIR3::wu];
    acc[dIR3::wv] += m_ratio_[i]*del_i[dIR3::wv];
    acc[dIR3::ww] += m_ratio_[i]*del_i[dIR3::ww];
  }
  return acc;
}
template<size_t N>
IR3 linear_combo_c1<N>::partial_t_contravariant(
    const IR3& position, double time) const {
  IR3 acc = {0, 0, 0};
  for (size_t i = 0; i < N; i++) acc +=
      m_ratio_[i]*field_set_[i]->partial_t_contravariant(
          position, t_ratio_[i]*time);
  return acc;
}

} // end namespace gyronimo.

#endif // GYRONIMO_LINEAR_COMBO_C1
