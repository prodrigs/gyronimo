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

// @morphism_cartesian.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_cartesian.hh>

namespace gyronimo {

morphism_cartesian::morphism_cartesian(double Lref) 
	: morphism(), Lref_(Lref), iLref_(1/Lref), Lref_3_(Lref*Lref*Lref) {
	
}
IR3 morphism_cartesian::operator()(const IR3 &q) const {
	return {Lref_*q[IR3::u], Lref_*q[IR3::v], Lref_*q[IR3::w]};
}
IR3 morphism_cartesian::inverse(const IR3 &x) const {
	return {iLref_*x[IR3::u], iLref_*x[IR3::v], iLref_*x[IR3::w]};
}
dIR3 morphism_cartesian::del(const IR3 &q) const {
	return {
		Lref_, 0.0, 0.0,
		0.0, Lref_, 0.0,
		0.0, 0.0, Lref_
	};
}
ddIR3 morphism_cartesian::ddel(const IR3 &q) const {
	return {
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	};
}
double morphism_cartesian::jacobian(const IR3 &q) const {
	return Lref_3_;
}
dIR3 morphism_cartesian::del_inverse(const IR3 &q) const {
	return {
		iLref_, 0.0, 0.0,
		0.0, iLref_, 0.0,
		0.0, 0.0, iLref_
	};
}
IR3 morphism_cartesian::to_covariant(const IR3 &q, const IR3 &A) const {
	return A;
}
IR3 morphism_cartesian::to_contravariant(const IR3 &q, const IR3 &A) const {
	return A;
}
IR3 morphism_cartesian::from_covariant(const IR3 &q, const IR3 &A) const {
	return A;
}
IR3 morphism_cartesian::from_contravariant(const IR3 &q, const IR3 &A) const {
	return A;
}

} // end namespace gyronimo