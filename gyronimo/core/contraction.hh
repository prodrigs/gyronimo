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

enum variance_type { covariant, contravariant };
enum contraction_index { first, second, third };

double inner_product(const IR3& A, const IR3& B);
IR3 contraction(const SM3& g, const IR3& B);
template<contraction_index> IR3 contraction(const dIR3& A, const IR3& B);
template<contraction_index> dIR3 contraction(const dSM3& A, const IR3& B);
template<contraction_index> dIR3 contraction(const dIR3& A, const SM3& B);
template<contraction_index> ddIR3 contraction(const ddIR3& A, const SM3& g);
template<contraction_index> ddIR3 contraction(const dIR3& A, const ddIR3& B);
dSM3 contraction(const SM3& g, const dSM3& d, const SM3& h);
IR3 cross_product(const IR3& A, const IR3& B);
template<variance_type>
IR3 cross_product(const IR3& A, const IR3& B, double jacobian);

}  // end namespace gyronimo.

#endif  // end GYRONIMO_CONTRACTION
