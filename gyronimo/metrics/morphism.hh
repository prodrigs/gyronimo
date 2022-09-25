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

// @morphism.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MORPHISM
#define GYRONIMO_MORPHISM

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/core/contraction.hh>

namespace gyronimo {

//! Abstract morphism class from curvilinear @f$ q^\alpha @f$ into cartesian @f$ \textbf{x} @f$ coordinates.
class morphism {

public:

	morphism() {};
	virtual ~morphism() {};

	//! Maps curvilinear coordinates @f$ q^\alpha @f$ into cartesian coordinates @f$ \textbf{x} @f$.
	virtual IR3 operator()(const IR3 &q) const = 0;
	//! Inverse transform from cartesian coordinates @f$ \textbf{x} @f$ into curvilinear coordinates @f$ q^\alpha @f$.
	virtual IR3 inverse(const IR3 &x) const = 0;
	//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point @f$ q^\alpha @f$.
	/*!
		Implements the covariant basis vectors:
		@f$ \textbf{e}_\alpha = \frac{\partial \textbf{x}}{\partial q^\alpha} @f$
	*/
	virtual dIR3 del(const IR3 &q) const = 0;
	//! Returns the morphism's second derivatives, calculated in point `q`.
	/*!
		Implements the derivatives of the covariant basis vectors:
		@f$ \partial_\beta \, \textbf{e}_\alpha = \frac{\partial^2 \textbf{x}}{\partial q^\beta \, \partial q^\alpha} @f$
	*/
	virtual ddIR3 del2(const IR3 &q) const = 0;

	virtual double jacobian(const IR3 &q) const;
	virtual dIR3 del_inverse(const IR3 &q) const;
	virtual dIR3 tan_basis(const IR3& q) const;
	virtual dIR3 dual_basis(const IR3& q) const;

	virtual IR3 to_covariant(const IR3 &q, const IR3 &A) const;
	virtual IR3 to_contravariant(const IR3 &q, const IR3 &A) const;
	virtual IR3 from_covariant(const IR3 &q, const IR3 &A) const;
	virtual IR3 from_contravariant(const IR3 &q, const IR3 &A) const;

private:

}; // end class morphism

} // end namespace gyronimo

#endif // GYRONIMO_MORPHISM
