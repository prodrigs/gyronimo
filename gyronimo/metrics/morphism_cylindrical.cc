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

// @morphism_cylindrical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_cylindrical.hh>

namespace gyronimo {

//! Maps the curvilinear coordinates `q = (r, \phi, z)` into cartesian coordinates `x`.
IR3 morphism_cylindrical::operator()(const IR3 &q) const {
	return {q[IR3::u] * std::cos(q[IR3::v]), q[IR3::u] * std::sin(q[IR3::v]), q[IR3::w]};
}

//! Inverse transform from cartesian coordinates `x` into curvilinear coordinates `q`.
IR3 morphism_cylindrical::inverse(const IR3 &x) const {
	return {sqrt(x[IR3::u]*x[IR3::u] + x[IR3::v]*x[IR3::v]),
			atan2(x[IR3::v], x[IR3::u]), x[IR3::w]};
}

//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point `q`.
dIR3 morphism_cylindrical::del(const IR3 &q) const {
	double sn = sin(q[IR3::v]);
	double cn = cos(q[IR3::v]);
	double r  = q[IR3::u];
	return {
		cn, - r * sn, 0,
		sn,   r * cn, 0,
		0, 		0, 	  1
	};
}

//! General-purpose implementation of the Jacobian of the transformation in point `q`.
double morphism_cylindrical::jacobian(const IR3 &q) const {
	return q[IR3::u];
}

//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
dIR3 morphism_cylindrical::del_inverse(const IR3 &q) const {
	double sn = sin(q[IR3::v]);
	double cn = cos(q[IR3::v]);
	double ir = 1 / q[IR3::u];
	return {
		cn, sn, 0,
		-sn * ir,  cn * ir, 0,
		0,	0   , 1
	};
}

} // end namespace gyronimo