#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

// Returns the jacobian of the transformation in point `q`.
double morphism::jacobian(const IR3 &q) const {
	dIR3 e = del(q);
	return (
		e[dIR3::uu] * (e[dIR3::vv] * e[dIR3::ww] - e[dIR3::vw] * e[dIR3::wv]) +
		e[dIR3::uv] * (e[dIR3::vw] * e[dIR3::wu] - e[dIR3::vu] * e[dIR3::ww]) +
		e[dIR3::uw] * (e[dIR3::vu] * e[dIR3::wv] - e[dIR3::vv] * e[dIR3::wu])
	);
}

// Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
dIR3 morphism::del_inverse(const IR3 &q) const {
	return gyronimo::inverse(del(q));
}

// Returns the covariant components of `A` in position `q`.
IR3 morphism::to_covariant(const IR3 &q, const IR3 &A) const {
	return contraction<first>(del(q), A);
}

// Returns the contravariant components of `A` in the position `q`.
IR3 morphism::to_contravariant(const IR3 &q, const IR3 &A) const {
	return contraction<second>(del_inverse(q), A);
}

// Returns the cartesian vector from its covariant components `A` in position `q`.
IR3 morphism::from_covariant(const IR3 &q, const IR3 &A) const {
	return contraction<first>(del_inverse(q), A);
}

// Returns the cartesian vector from its contravariant components `A` in position `q`.
IR3 morphism::from_contravariant(const IR3 &q, const IR3 &A) const {
	return contraction<second>(del(q), A);
}

} // end namespace gyronimo