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

// @equilibrium_circular.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_EQUILIBRIUM_CIRCULAR
#define GYRONIMO_EQUILIBRIUM_CIRCULAR

#include <gyronimo/fields/IR3field_c1.hh>
#include <gyronimo/metrics/metric_polar_torus.hh>

namespace gyronimo {

//! Static toroidal equilibrium with centred circular magnetic surfaces.
/*!
    Toroidal coordinates are defined by the `metric_polar_torus` object and the
    poloidal flux is constant over circumferences of constant `r`. Poloidal
    dependence comes from the 1/R toroidal term only. The radial dependence is
    established by the safety-factor (`q` and `qprime`), supplied as lambda
    objects of the type `radial_profile`. Inherited member functions
    `covariant`, `magnitude`, `covariant_versor`, `contravariant_versor`,
    `curl`, `del_magnitude` and `partial_t_magnitude` are reimplemented for
    efficiency purposes.
*/
class equilibrium_circular : public IR3field_c1 {
 public:
  typedef std::function<double(double)> radial_profile;

  equilibrium_circular(
      double m_factor, const metric_polar_torus *g,
      radial_profile q, radial_profile qprime)
      : IR3field_c1(m_factor, 1.0, g), metric_(g),
        q_(std::move(q)), qprime_(std::move(qprime)) {};
  virtual ~equilibrium_circular() override {};

  virtual IR3 contravariant(const IR3& position, double time) const override;
  virtual dIR3 del_contravariant(
      const IR3& position, double time) const override;
  virtual IR3 partial_t_contravariant(
      const IR3& position, double time) const override {return {0.0,0.0,0.0};};

  virtual IR3 covariant(const IR3& position, double time) const override;
  virtual double magnitude(const IR3& position, double time) const override;
  virtual IR3 covariant_versor(const IR3& position, double time) const override;
  virtual IR3 contravariant_versor(
      const IR3& position, double time) const override;

  virtual IR3 del_magnitude(const IR3& position, double time) const override;
  virtual double partial_t_magnitude(
      const IR3& position, double time) const override {return 0.0;};
  virtual IR3 curl(const IR3& position, double time) const override;

  double q(double r) const {return q_(r);};
  double qprime(double r) const {return qprime_(r);};

 private:
  const metric_polar_torus *metric_;
  const radial_profile q_;
  const radial_profile qprime_;
  };

} // end namespace gyronimo.

#endif // GYRONIMO_EQUILIBRIUM_CIRCULAR
