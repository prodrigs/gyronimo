// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Manuel Assunção.

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

// @metric_cylindrical.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_CYLINDRICAL
#define GYRONIMO_METRIC_CYLINDRICAL

#include <gyronimo/metrics/metric_nexus.hh>
#include <gyronimo/metrics/morphism_cylindrical.hh>

namespace gyronimo {

//! Trivial covariant metric for cylindrical space.
class metric_cylindrical : public metric_nexus {
public:
	metric_cylindrical(const morphism_cylindrical *morph);
	virtual ~metric_cylindrical() override {};

	virtual SM3 operator()(const IR3& q) const override;
	virtual SM3 inverse(const IR3& q) const override;

	virtual dSM3 del(const IR3& q) const override;
	virtual double jacobian(const IR3& q) const override;
	virtual IR3 del_jacobian(const IR3& q) const override;
	virtual IR3 to_covariant(const IR3& B, const IR3& q) const override;
	virtual IR3 to_contravariant(const IR3& B, const IR3& q) const override;

	virtual ddIR3 christoffel_first_kind(const IR3& q) const override;
	virtual ddIR3 christoffel_second_kind(const IR3& q) const override;

	double Lref() {return Lref_;};

private:
	const double Lref_, Lref_2_, iLref_2_, Lref_3_;

}; // end class metric_cylindrical

} // end namespace gyronimo.

#endif // GYRONIMO_METRIC_CARTESIAN
