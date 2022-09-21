#include <gyronimo/metrics/morphism_spherical.hh>

namespace gyronimo {

// Maps the curvilinear coordinates `q` into cartesian coordinates `x`.
IR3 morphism_spherical::operator()(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::w]);
	double sn_theta = std::sin(q[IR3::w]);
	double cn_phi = std::cos(q[IR3::v]);
	double sn_phi = std::sin(q[IR3::v]);
	return {r * cn_phi * sn_theta, r * sn_phi * sn_theta, r * cn_theta};
}

// Inverse transform from cartesian coordinates `x` into curvilinear coordinates `q`.
IR3 morphism_spherical::inverse(const IR3 &x) const {
	double rho = x[IR3::u] * x[IR3::u] + x[IR3::v] * x[IR3::v];
	return {std::sqrt(rho + x[IR3::w] * x[IR3::w]), std::atan2(x[IR3::u], x[IR3::v]), std::atan2(x[IR3::w], sqrt(rho))};
}

// Returns the morphism's derivatives, correspondent to the covariant basis vectors in point `q`.
dIR3 morphism_spherical::del(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::w]);
	double sn_theta = std::sin(q[IR3::w]);
	double cn_phi = std::cos(q[IR3::v]);
	double sn_phi = std::sin(q[IR3::v]);
	return {
		cn_phi * sn_theta, - r * sn_phi * sn_theta, r * cn_phi * cn_theta,
		sn_phi * sn_theta,   r * cn_phi * sn_theta, r * sn_phi * cn_theta,
				 cn_theta, 						 0, - r * sn_theta
	};
}

dSM3 morphism_spherical::g_del(const IR3 &q) const {
	double r = q[IR3::u];
	double sn_theta = std::sin(q[IR3::w]);
	double sn_2theta = std::sin(2 * q[IR3::w]);
	return {
		0, 0, 0, // g_uu
		0, 0, 0, // g_uv
		0, 0, 0, // g_uw
		2 * r * sn_theta * sn_theta, 0, r * r * sn_2theta, // g_vv
		0, 0, 0, // g_vw
		2 * r, 0, 0  // g_ww
	};
}

// Returns the jacobian of the transformation in point `q`.
double morphism_spherical::jacobian(const IR3 &q) const {
	return q[IR3::u] * q[IR3::u] * sin(IR3::w);
}

// Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
dIR3 morphism_spherical::del_inverse(const IR3 &q) const {
	double ir = 1 / q[IR3::u];
	double cn_theta = std::cos(q[IR3::w]);
	double sn_theta = std::sin(q[IR3::w]);
	double csc_theta = 1 / sn_theta;
	double cn_phi = std::cos(q[IR3::v]);
	double sn_phi = std::sin(q[IR3::v]);
	return {
			   cn_phi * sn_theta, 		sn_phi * sn_theta, 		cn_theta,
		- ir * sn_phi * csc_theta, ir * cn_phi * csc_theta, 			0,
		  ir * cn_phi * cn_theta,  ir * sn_phi * cn_theta, - ir * sn_theta
	};
}

} // end namespace gyronimo