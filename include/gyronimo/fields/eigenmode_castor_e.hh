// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @eigenmode_castor_e.hh

#ifndef GYRONIMO_EIGENMODE_CASTOR_E
#define GYRONIMO_EIGENMODE_CASTOR_E

#include <gyronimo/fields/eigenmode_castor_a.hh>

namespace gyronimo {

//! Electric-field eigenvector from `CASTOR` ceig output file.
/*!
    The electric-field **covariant** components are extracted as @f$\mathbf{E} =
    - \partial_t \mathbf{A} = i \omega \mathbf{A} @f$ from the vector potential
    stored in a `eigenmode_castor_a`, with a similar normalisation to the
    on-axis Alfven time. The normalisation of the electric-field magnitude is
    set by `m_factor` (SI).
*/
class eigenmode_castor_e : public IR3field {
 public:
  eigenmode_castor_e(
      double m_factor, double t_factor,
      const parser_castor *p, const metric_helena *g,
      const interpolator1d_factory* ifactory)
      : IR3field(m_factor, t_factor, g),
        iw_(-p->w_imag(), p->w_real()),
        A_(m_factor*t_factor, t_factor, p, g, ifactory) {};
  virtual ~eigenmode_castor_e() override {};

  virtual IR3 covariant(const IR3& position, double time) const override;
  virtual IR3 contravariant(const IR3& position, double time) const override;

  const parser_castor* parser() const {return A_.parser();};

 private:
  std::complex<double> iw_;
  eigenmode_castor_a A_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_EIGENMODE_CASTOR_E
