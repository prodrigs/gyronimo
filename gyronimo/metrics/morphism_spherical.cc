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

//! Maps spherical coordinates @f$ q^\alpha @f$ into cartesian coordinates @f$ \textbf{x} @f$.
/*!
	In spherical coordinates, @f$ \left( q^1, q^2, q^3 \right) = \left( r, \theta, \phi \right) @f$.
	Implements the coordinate transformation:
	@f{gather*}{
		\begin{aligned}
			x &= r \, \cos \phi \, \sin \theta \\
			y &= r \, \sin \phi \, \sin \theta \\
			z &= r \, \cos \theta
		\end{aligned}
    @f}
*/
IR3 morphism_spherical::operator()(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::v]);
	double sn_theta = std::sin(q[IR3::v]);
	double cn_phi = std::cos(q[IR3::w]);
	double sn_phi = std::sin(q[IR3::w]);
	return {r * cn_phi * sn_theta, r * sn_phi * sn_theta, r * cn_theta};
}

//! Inverse transform from cartesian coordinates @f$ \textbf{x} @f$ into spherical coordinates @f$ q^\alpha @f$.
/*!
	Implements the inverse transformation:
	@f{gather*}{
		\begin{aligned}
			r &= \sqrt{x^2 + y^2 + z^2} \\
			\theta &= \arctan \left( \frac{\sqrt{x^2 + y^2}}{z} \right) \\
			\phi &= \arctan \left( \frac{y}{x} \right)
		\end{aligned}
    @f}
*/
IR3 morphism_spherical::inverse(const IR3 &x) const {
	double rho = x[IR3::u] * x[IR3::u] + x[IR3::v] * x[IR3::v];
	return {std::sqrt(rho + x[IR3::w] * x[IR3::w]), std::atan2(sqrt(rho), x[IR3::w]), std::atan2(x[IR3::v], x[IR3::u])};
}

//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point @f$ q^\alpha @f$.
/*!
	Implements the coordinate transformation's first derivatives:
	@f{gather*}{
		\begin{aligned}
			\textbf{e}_r &= \left( \cos \phi \, \sin \theta, \sin \phi \, \sin \theta, \cos \theta \right) \\
			\textbf{e}_\theta &= \left( r \, \cos \phi \, \cos \theta, r \, \sin \phi \, \cos \theta, -r \, \sin \theta \right) \\
			\textbf{e}_\phi &= \left( -r \, \sin \phi \, \sin \theta, r\, \cos \phi \, \sin \theta, 0 \right)
		\end{aligned}
    @f}
*/
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

//! Returns the morphism's second derivatives, calculated in point @f$ q^\alpha @f$.
ddIR3 morphism_spherical::del2(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::w]);
	double sn_theta = std::sin(q[IR3::w]);
	double cn_phi = std::cos(q[IR3::v]);
	double sn_phi = std::sin(q[IR3::v]);
	return {
	// iuu		iuv					iuw						ivv						ivw						iww
		0, cn_phi * cn_theta, -sn_phi * sn_theta, -r * cn_phi * sn_theta, -r * sn_phi * cn_theta, -r * cn_phi * sn_theta,
		0, sn_phi * cn_theta,  cn_phi * sn_theta, -r * sn_phi * sn_theta,  r * cn_phi * cn_theta, -r * sn_phi * sn_theta,
		0, 		   -sn_theta,				   0, 		   -r * cn_theta,					   0,		0
	};
}

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
/*!
	Implements the Jacobian in spherical coordinates: 
	@f$ J = r^2 \, \sin \theta @f$
*/
double morphism_spherical::jacobian(const IR3 &q) const {
	return q[IR3::u] * q[IR3::u] * std::sin(q[IR3::v]);
}

//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point @f$ q^\alpha @f$.
/*!
	Implements the inverse transformation's first derivatives:
	@f{gather*}{
		\begin{aligned}
			\textbf{e}^r &= \left( \cos \phi \, \sin \theta, \sin \phi \, \sin \theta, \cos \theta \right) \\
			\textbf{e}^\theta &= \left( \frac{1}{r} \, \cos \phi \, \cos \theta, \frac{1}{r} \, \sin \phi \, \cos \theta, -\frac{1}{r} \, \sin \theta \right) \\
			\textbf{e}^\phi &= \left( -\frac{\sin \phi}{r \, \sin \theta}, \frac{\cos \phi}{r \, \sin \theta}, 0 \right)
		\end{aligned}
    @f}
*/
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