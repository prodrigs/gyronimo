// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues and Manuel Assunção.

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

// @morphism.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
/*!
	Implements the rule:
	@f$ J = \textbf{e}_1 \cdot \left( \textbf{e}_2 \times \textbf{e}_3 \right) @f$
*/
double morphism::jacobian(const IR3 &q) const {
	dIR3 e = del(q);
	return (
		e[dIR3::uu] * (e[dIR3::vv] * e[dIR3::ww] - e[dIR3::vw] * e[dIR3::wv]) +
		e[dIR3::uv] * (e[dIR3::vw] * e[dIR3::wu] - e[dIR3::vu] * e[dIR3::ww]) +
		e[dIR3::uw] * (e[dIR3::vu] * e[dIR3::wv] - e[dIR3::vv] * e[dIR3::wu])
	);
}

//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point @f$ q^\alpha @f$.
/*!
	Implements the covariant basis vectors:
	@f$ \textbf{e}^\alpha = \nabla q^\alpha @f$
*/
dIR3 morphism::del_inverse(const IR3 &q) const {
	return gyronimo::inverse(del(q));
}

//! Returns the tangent basis @f$ \textbf{e}_\alpha @f$
dIR3 morphism::tan_basis(const IR3& q) const {
	return del(q);
}

//! Returns the dual basis @f$ \textbf{e}^\alpha @f$
dIR3 morphism::dual_basis(const IR3& q) const {
	return del_inverse(q);
}

//! Returns the covariant components of @f$ \textbf{A} @f$ in position @f$ q^\alpha @f$.
/*!
	Implements the rule:
	@f$ A_\alpha = \textbf{A} \cdot \textbf{e}_\alpha @f$
*/
IR3 morphism::to_covariant(const IR3 &A, const IR3 &q) const {
	return contraction<first>(del(q), A);
}

//! Returns the contravariant components of @f$ \textbf{A} @f$ in the position @f$ q^\alpha @f$.
/*!
	Implements the rule:
	@f$ A^\alpha = \textbf{A} \cdot \textbf{e}^\alpha @f$
*/
IR3 morphism::to_contravariant(const IR3 &A, const IR3 &q) const {
	return contraction<second>(del_inverse(q), A);
}

//! Returns the cartesian vector from its covariant components @f$ A_\alpha @f$ in position @f$ q^\alpha @f$.
/*!
	Implements the rule:
	@f$ \textbf{A} = A_\alpha \, \textbf{e}^\alpha @f$
*/
IR3 morphism::from_covariant(const IR3 &A, const IR3 &q) const {
	return contraction<first>(del_inverse(q), A);
}

//! Returns the cartesian vector from its contravariant components @f$ A^\alpha @f$ in position @f$ q^\alpha @f$.
/*!
	Implements the rule:
	@f$ \textbf{A} = A^\alpha \, \textbf{e}_\alpha @f$
*/
IR3 morphism::from_contravariant(const IR3 &A, const IR3 &q) const {
	return contraction<second>(del(q), A);
}

} // end namespace gyronimo