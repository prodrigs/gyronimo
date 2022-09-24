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

//! Abstract morphism from curvilinear (q) into cartesian (x) coordinates.
class morphism {

public:

	morphism() {};
	virtual ~morphism() {};

	//! Maps the curvilinear coordinates `q` into cartesian coordinates `x`.
	virtual IR3 operator()(const IR3 &q) const = 0;
	//! Inverse transform from cartesian coordinates `x` into curvilinear coordinates `q`.
	virtual IR3 inverse(const IR3 &x) const = 0;
	//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point `q`.
	virtual dIR3 del(const IR3 &q) const = 0;

	virtual double jacobian(const IR3 &q) const;
	virtual dIR3 del_inverse(const IR3 &q) const;
	virtual dIR3 tan_basis(const IR3& q) const {return del(q);};
	virtual dIR3 dual_basis(const IR3& q) const {return del_inverse(q);};

	virtual IR3 to_covariant(const IR3 &q, const IR3 &A) const;
	virtual IR3 to_contravariant(const IR3 &q, const IR3 &A) const;
	virtual IR3 from_covariant(const IR3 &q, const IR3 &A) const;
	virtual IR3 from_contravariant(const IR3 &q, const IR3 &A) const;

private:

}; // end class morphism

} // end namespace gyronimo

#endif // GYRONIMO_MORPHISM
