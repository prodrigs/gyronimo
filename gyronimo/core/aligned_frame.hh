// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2024 Paulo Rodrigues.

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

// @aligned_frame.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_ALIGNED_FRAME
#define GYRONIMO_ALIGNED_FRAME

#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/metric_connected.hh>

namespace gyronimo {

//! Defines a reference frame aligned with a given `IR3field` at every point.
/*!
    At a given position `q` where the underlying `IR3field` evaluates to `B(q)`,
    the aligned frame is produced as follows:

    1. Rotate the cartesian frame {ux, uy, uz}->{ux', uy', uz} around the uz
    axis so that uy' lies in the plane defined by uz and `B(q)`;

    2. Rotate again {ux', uy', uz}->{ux', uy'', uz'} around the ux' axis in
    order to align uz' with `B(q)`;

    Besides its versors and gyro_versors [i.e., the vectors @f$\{\mathbf{a} =
    \mathbf{b} \times \mathbf{c}@f$ as defined in Littlejohn, J. Plasma Phys.
    **29**, 111 (1983)], this class provides also functions to convert a
    cartesian velocity (assumed to be normalised to some v_ref) to a set {energy
    (normalised to @f$\frac{1}{2} m v_\textrm{ref}^2@f$), pitch} and back.
    Whenever needed, the gyrophase is defined as the angle around the magnetic
    field (or uz' axis) measured in the ux'-uy'' plane with origin at ux'.
*/
class aligned_frame {
 public:
  struct energy_data_t {
    double energy, pitch, charge_sign;
  };
  aligned_frame(const IR3field* field);

  std::array<IR3, 3> versors(const IR3& q, double t) const;
  std::array<IR3, 3> gyro_versors(
      double gyrophase, const IR3& q, double t) const;

  IR3 velocity_from_energy_data(
      const energy_data_t& data, double gyrophase, const IR3& q,
      double t) const;
  energy_data_t energy_data_from_velocity(
      const IR3& v, double charge_sign, const IR3& q, double t) const;
 private:
  const IR3field* field_;
  const metric_connected* metric_;
};

}  // end namespace gyronimo.

#endif  // GYRONIMO_ALIGNED_FRAME
