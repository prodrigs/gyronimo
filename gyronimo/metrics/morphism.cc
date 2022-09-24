// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

//! General-purpose implementation of the Jacobian of the transformation in point `q`.
double morphism::jacobian(const IR3 &q) const {
	dIR3 e = del(q);
	return (
		e[dIR3::uu] * (e[dIR3::vv] * e[dIR3::ww] - e[dIR3::vw] * e[dIR3::wv]) +
		e[dIR3::uv] * (e[dIR3::vw] * e[dIR3::wu] - e[dIR3::vu] * e[dIR3::ww]) +
		e[dIR3::uw] * (e[dIR3::vu] * e[dIR3::wv] - e[dIR3::vv] * e[dIR3::wu])
	);
}

//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
dIR3 morphism::del_inverse(const IR3 &q) const {
	return gyronimo::inverse(del(q));
}

//! Returns the covariant components of `A` in position `q`.
IR3 morphism::to_covariant(const IR3 &q, const IR3 &A) const {
	return contraction<first>(del(q), A);
}

//! Returns the contravariant components of `A` in the position `q`.
IR3 morphism::to_contravariant(const IR3 &q, const IR3 &A) const {
	return contraction<second>(del_inverse(q), A);
}

//! Returns the cartesian vector from its covariant components `A` in position `q`.
IR3 morphism::from_covariant(const IR3 &q, const IR3 &A) const {
	return contraction<first>(del_inverse(q), A);
}

//! Returns the cartesian vector from its contravariant components `A` in position `q`.
IR3 morphism::from_contravariant(const IR3 &q, const IR3 &A) const {
	return contraction<second>(del(q), A);
}

} // end namespace gyronimo