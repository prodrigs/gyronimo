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

// @morphism_spherical.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_SPHERICAL
#define GYRONIMO_MORPHISM_SPHERICAL

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

class morphism_spherical : public morphism {

public:

	morphism_spherical(const double &r0) 
		: morphism(), r0_(r0), ir0_(1/r0), r0_3_(r0*r0*r0) {};
	virtual ~morphism_spherical() override {};

	virtual IR3 operator()(const IR3 &q) const override;
	virtual IR3 inverse(const IR3 &x) const override;
	virtual dIR3 del(const IR3 &q) const override;
	virtual ddIR3 ddel(const IR3 &q) const override;

	virtual double jacobian(const IR3 &q) const override;
	virtual dIR3 del_inverse(const IR3 &q) const override;

	double r0() {return r0_;};

private:
	const double r0_, ir0_, r0_3_;

}; // end class morphism_spherical

} // end namespace gyronimo

#endif // GYRONIMO_MORPHISM_SPHERICAL