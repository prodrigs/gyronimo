// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues.

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

#include <cmath>

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

class morphism_spherical : public morphism {

public:

	morphism_spherical() : morphism() {};
	virtual ~morphism_spherical() override {};

	virtual IR3 operator()(const IR3 &q) const override;
	virtual IR3 inverse(const IR3 &x) const override;
	virtual dIR3 del(const IR3 &q) const override;
	virtual ddIR3 del2(const IR3 &q) const override;

	virtual double jacobian(const IR3 &q) const override;
	virtual dIR3 del_inverse(const IR3 &q) const override;

private:

}; // end class morphism_spherical

} // end namespace gyronimo

#endif // GYRONIMO_MORPHISM_SPHERICAL