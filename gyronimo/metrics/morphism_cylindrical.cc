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

// @morphism_cylindrical.cc, this file is part of ::gyronimo::

#include <cmath>
#include <gyronimo/metrics/morphism_cylindrical.hh>

namespace gyronimo {

//! Maps cylindrical coordinates @f$ q^\alpha @f$ into cartesian coordinates @f$ \textbf{x} @f$.
/*!
	In cylindrical coordinates, @f$ \left( q^1, q^2, q^3 \right) = \left( R, \phi, Z \right) @f$.
	Implements the coordinate transformation:
	@f{gather*}{
		\begin{aligned}
			x &= L_{ref} \, R \, \cos \phi \\
			y &= L_{ref} \, R \, \sin \phi \\
			z &= L_{ref} \, Z
		\end{aligned}
    @f}
*/
IR3 morphism_cylindrical::operator()(const IR3 &q) const {
	return {Lref_ * q[IR3::u] * std::cos(q[IR3::v]), 
			Lref_ * q[IR3::u] * std::sin(q[IR3::v]), 
			Lref_ * q[IR3::w]};
}

//! Inverse transform from cartesian coordinates @f$ \textbf{x} @f$ into cylindrical coordinates @f$ q^\alpha @f$.
/*!
	Implements the inverse transformation:
	@f{gather*}{
		\begin{aligned}
			R &= \frac{1}{L_{ref}} \sqrt{x^2 + y^2} \\
			\phi &= \arctan \left( \frac{y}{x} \right) \\
			Z &= \frac{z}{L_{ref}}
		\end{aligned}
    @f}
*/
IR3 morphism_cylindrical::inverse(const IR3 &x) const {
	return {iLref_ * std::sqrt(x[IR3::u]*x[IR3::u] + x[IR3::v]*x[IR3::v]),
			std::atan2(x[IR3::v], x[IR3::u]), 
			iLref_ * x[IR3::w]};
}

//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point @f$ q^\alpha @f$.
/*!
	Implements the coordinate transformation's first derivatives:
	@f{gather*}{
		\begin{aligned}
			\textbf{e}_R &= \left( L_{ref} \, \cos \phi, L_{ref} \, \sin \phi, 0 \right) \\
			\textbf{e}_\phi &= \left( -L_{ref} \, R \, \sin \phi, L_{ref} \, R \, \cos \phi, 0 \right) \\
			\textbf{e}_Z &= \left( 0, 0, L_{ref} \right)
		\end{aligned}
    @f}
*/
dIR3 morphism_cylindrical::del(const IR3 &q) const {
	double sn = Lref_ * std::sin(q[IR3::v]);
	double cn = Lref_ * std::cos(q[IR3::v]);
	double r  = q[IR3::u];
	return {
		cn, - r * sn, 0,
		sn,   r * cn, 0,
		0, 		0, 	 Lref_
	};
}

//! Returns the morphism's second derivatives, calculated in point @f$ q^\alpha @f$.
ddIR3 morphism_cylindrical::ddel(const IR3 &q) const {
	double r = q[IR3::u];
	double sn = Lref_ * std::sin(q[IR3::v]);
	double cn = Lref_ * std::cos(q[IR3::v]);
	return {
	// iuu iuv iuw   ivv  ivw  iww
		0, -sn, 0, -r * cn, 0, 0,
		0,  cn, 0, -r * sn, 0, 0,
		0,   0, 0,		 0, 0, 0
	};
}

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
/*!
	Implements the Jacobian in cylindrical coordinates: 
	@f$ J = {L_{ref}}^3 \, R @f$
*/
double morphism_cylindrical::jacobian(const IR3 &q) const {
	return Lref_3_ * q[IR3::u];
}

//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point @f$ q^\alpha @f$.
/*!
	Implements the inverse transformation's first derivatives:
	@f{gather*}{
		\begin{aligned}
			\textbf{e}^R &= \left( \frac{1}{L_{ref}} \, \cos \phi, \frac{1}{L_{ref}} \, \sin \phi, 0 \right) \\
			\textbf{e}^\phi &= \left( -\frac{1}{L_{ref} \, r} \, \sin \phi, \frac{1}{L_{ref} \, r} \, \cos \phi, 0 \right) \\
			\textbf{e}^Z &= \left( 0, 0, \frac{1}{L_{ref}} \right)
		\end{aligned}
    @f}
*/
dIR3 morphism_cylindrical::del_inverse(const IR3 &q) const {
	double sn = iLref_ * std::sin(q[IR3::v]);
	double cn = iLref_ * std::cos(q[IR3::v]);
	double ir = 1 / q[IR3::u];
	return {
		cn, sn, 0,
		-sn * ir,  cn * ir, 0,
		0,	0, iLref_
	};
}

} // end namespace gyronimo