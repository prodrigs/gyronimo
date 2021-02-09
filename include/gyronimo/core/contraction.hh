// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @contraction.hh

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
template<contraction_index>
dIR3 contraction(const dSM3& dA, const IR3& B);
template<contraction_index, contraction_index>
dIR3 contraction(const SM3& g, const dIR3& dB);

} // end namespace gyronimo.

#endif // end GYRONIMO_CONTRACTION
