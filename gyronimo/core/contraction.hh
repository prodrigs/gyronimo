// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2023 Paulo Rodrigues and Manuel Assunção.

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

// @contraction.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CONTRACTION
#define GYRONIMO_CONTRACTION

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/core/SM3algebra.hh>

namespace gyronimo {

enum variance_type {covariant, contravariant};
enum contraction_index {first, second, third};

IR3 contraction(const SM3& g, const IR3& B);
double inner_product(const IR3& A, const IR3& B);
IR3 cross_product(const IR3& A, const IR3& B);
dSM3 contraction(const SM3& g, const dSM3& d, const SM3& h);

template<variance_type>
IR3 cross_product(const IR3& A, const IR3& B, double jacobian);
template<contraction_index>
IR3 contraction(const dIR3& dA, const IR3& B);
template<contraction_index>
dIR3 contraction(const dSM3& dA, const IR3& B);
template<contraction_index, contraction_index>
IR3 contraction(const ddIR3& ddA, const IR3& B, const IR3& C);
template<contraction_index, contraction_index>
ddIR3 contraction(const SM3& g, const ddIR3& ddA);
template<contraction_index, contraction_index>
ddIR3 contraction(const dIR3& dA, const ddIR3& ddB);
template<contraction_index, contraction_index>
dIR3 contraction(const SM3& g, const dIR3& dB);

} // end namespace gyronimo.

#endif // end GYRONIMO_CONTRACTION
