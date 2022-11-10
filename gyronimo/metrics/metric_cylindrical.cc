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

// @metric_cylindrical.cc, this file is part of ::gyronimo::

#include <cmath>
#include <gyronimo/metrics/metric_cylindrical.hh>

namespace gyronimo {

metric_cylindrical::metric_cylindrical(const morphism_cylindrical *morph) 
		: metric_connected(morph),
		Lref_(morph->Lref()), Lref_2_(Lref_*Lref_), 
		iLref_2_(1/Lref_2_), Lref_3_(Lref_*Lref_*Lref_) {
	
}

//! General-purpose implementation of the covariant metric.
/*!
	Implements the covariant metric in cylindrical coordinates: 
	@f$ g_{\alpha\beta} = \left(\begin{matrix}
		L_{ref}^2 && 0 && 0 \\
		0 && L_{ref}^2 \, R^2 && 0 \\
		0 && 0 && L_{ref}^2
	\end{matrix}\right) @f$
*/
SM3 metric_cylindrical::operator()(const IR3& q) const {
	return {Lref_2_, 0.0, 0.0, Lref_2_*q[IR3::u]*q[IR3::u], 0.0, Lref_2_};
}

//! General-purpose implementation of the inverse (i.e., contravariant metric).
/*!
	Implements the contravariant metric in cylindrical coordinates: 
	@f$ g^{\alpha\beta} = \left(\begin{matrix}
		\frac{1}{L_{ref}^2} && 0 && 0 \\
		0 && \frac{1}{L_{ref}^2 \, R^2} && 0 \\
		0 && 0 && \frac{1}{L_{ref}^2}
	\end{matrix}\right) @f$
*/
SM3 metric_cylindrical::inverse(const IR3& q) const {
	return {iLref_2_, 0.0, 0.0, iLref_2_/(q[IR3::u]*q[IR3::u]), 0.0, iLref_2_};
}

dSM3 metric_cylindrical::del(const IR3& q) const {
	return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 2*q[IR3::u]*Lref_2_, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
/*!
	Implements the Jacobian in cylindrical coordinates: 
	@f$ J = L_{ref}^3 \, R @f$
*/
double metric_cylindrical::jacobian(const IR3& q) const {
	return Lref_3_*q[IR3::u];
}

//! General-purpose implementation of the Jacobian gradient in point @f$ q^\alpha @f$.
/*!
	Implements the Jacobian gradient in cylindrical coordinates: 
	@f$ \nabla J = \left( L_{ref}^2\,\cos\phi, L_{ref}^2\,\sin\phi, 0 \right) @f$
*/
IR3 metric_cylindrical::del_jacobian(const IR3& q) const {
	return {Lref_3_, 0.0, 0.0};
}

IR3 metric_cylindrical::to_covariant(const IR3& B, const IR3& q) const {
	return {Lref_2_*B[IR3::u], Lref_2_*B[IR3::v]*q[IR3::u]*q[IR3::u], Lref_2_*B[IR3::w]};
}

IR3 metric_cylindrical::to_contravariant(const IR3& B, const IR3& q) const {
	return {iLref_2_*B[IR3::u], iLref_2_*B[IR3::v]/(q[IR3::u]*q[IR3::u]), iLref_2_*B[IR3::w]};
}

ddIR3 metric_cylindrical::christoffel_first_kind(const IR3& q) const {
	return {
		0, 0, 0, -Lref_2_*q[IR3::u], 0, 0,
		0,  Lref_2_*q[IR3::u], 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0
	};
}

ddIR3 metric_cylindrical::christoffel_second_kind(const IR3& q) const {
	return {
		0, 0, 0,  -q[IR3::u], 0, 0,
		0, 1/q[IR3::u], 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0
	};
}

} // end namespace gyronimo