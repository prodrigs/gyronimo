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

// @morphism_spherical.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_spherical.hh>

namespace gyronimo {

//! Maps the spherical coordinates `q = (r, \theta, \phi)` into cartesian coordinates `x`.
IR3 morphism_spherical::operator()(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::v]);
	double sn_theta = std::sin(q[IR3::v]);
	double cn_phi = std::cos(q[IR3::w]);
	double sn_phi = std::sin(q[IR3::w]);
	return {r * cn_phi * sn_theta, r * sn_phi * sn_theta, r * cn_theta};
}

//! Inverse transform from cartesian coordinates `x` into curvilinear coordinates `(r, \theta, \phi)`.
IR3 morphism_spherical::inverse(const IR3 &x) const {
	double rho = x[IR3::u] * x[IR3::u] + x[IR3::v] * x[IR3::v];
	return {std::sqrt(rho + x[IR3::w] * x[IR3::w]), std::atan2(sqrt(rho), x[IR3::w]), std::atan2(x[IR3::v], x[IR3::u])};
}

//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point `q`.
dIR3 morphism_spherical::del(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::w]);
	double sn_theta = std::sin(q[IR3::w]);
	double cn_phi = std::cos(q[IR3::v]);
	double sn_phi = std::sin(q[IR3::v]);
	return {
		cn_phi * sn_theta, r * cn_phi * cn_theta, - r * sn_phi * sn_theta,
		sn_phi * sn_theta, r * sn_phi * cn_theta,   r * cn_phi * sn_theta,
				 cn_theta, 		  - r * sn_theta, 	0
	};
}

//! General-purpose implementation of the Jacobian of the transformation in point `q`.
double morphism_spherical::jacobian(const IR3 &q) const {
	return q[IR3::u] * q[IR3::u] * std::sin(q[IR3::v]);
}

//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
dIR3 morphism_spherical::del_inverse(const IR3 &q) const {
	double ir = 1 / q[IR3::u];
	double cn_theta = std::cos(q[IR3::w]);
	double sn_theta = std::sin(q[IR3::w]);
	double csc_theta = 1 / sn_theta;
	double cn_phi = std::cos(q[IR3::v]);
	double sn_phi = std::sin(q[IR3::v]);
	return {
			   cn_phi * sn_theta, 		sn_phi * sn_theta, 		  cn_theta,
		  ir * cn_phi * cn_theta,  ir * sn_phi * cn_theta, - ir * sn_theta,
		- ir * sn_phi * csc_theta, ir * cn_phi * csc_theta, 	0
	};
}

} // end namespace gyronimo