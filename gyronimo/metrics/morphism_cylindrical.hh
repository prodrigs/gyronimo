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

// @morphism_cylindrical.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM_CYLINDRICAL
#define GYRONIMO_MORPHISM_CYLINDRICAL

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

class morphism_cylindrical : public morphism {

public:

	morphism_cylindrical(const double &L0) 
		: morphism(), L0_(L0), iL0_(1/L0), L0_3_(L0*L0*L0) {};
	virtual ~morphism_cylindrical() override {};

	virtual IR3 operator()(const IR3 &q) const override;
	virtual IR3 inverse(const IR3 &x) const override;
	virtual dIR3 del(const IR3 &q) const override;
	virtual ddIR3 ddel(const IR3 &q) const override;

	virtual double jacobian(const IR3 &q) const override;
	virtual dIR3 del_inverse(const IR3 &q) const override;

	double L0() const {return L0_;};

private:
	const double L0_, iL0_, L0_3_;

}; // end class morphism_cylindrical

} // end namespace gyronimo

#endif // GYRONIMO_MORPHISM_CYLINDRICAL