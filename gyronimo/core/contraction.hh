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

// @contraction.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CONTRACTION
#define GYRONIMO_CONTRACTION

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/core/SM3algebra.hh>

namespace gyronimo {

IR3 contraction(const SM3& g, const IR3& B);
IR3 cross_product(const IR3& A, const IR3& B, double jacobian);
double inner_product(const IR3& A, const IR3& B);

//! Index to contract in a templated `contraction` of multi-indexed objects.
enum contraction_index {first, second, third};

template<contraction_index>
IR3 contraction(const dIR3& dA, const IR3& B);
template<contraction_index, contraction_index>
IR3 contraction(const ddIR3& ddA, const IR3& B, const IR3& C);
template<contraction_index>
dIR3 contraction(const dSM3& dA, const IR3& B);
template<contraction_index, contraction_index>
dIR3 contraction(const SM3& g, const dIR3& dB);
dSM3 contraction(const SM3& g, const dSM3& d, const SM3& h);

} // end namespace gyronimo.

#endif // end GYRONIMO_CONTRACTION
