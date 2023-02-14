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

// @mrophism_vmec.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_vmec.hh>
#include <numbers>
#include <gyronimo/core/multiroot.hh>

namespace gyronimo{

morphism_vmec::morphism_vmec(const parser_vmec *p, 
		const interpolator1d_factory *ifactory)
		: parser_(p), xm_(p->xm()), xn_(p->xn()) {

	//set radial grid block
    dblock_adapter s_range(p->radius());
    // set spectral components 
    Rmnc_ = new interpolator1d* [xm_.size()];
    Zmns_ = new interpolator1d* [xm_.size()];

	for(size_t i = 0; i < xm_.size(); ++i) {
		std::slice s_cut (i, s_range.size(), xm_.size());
		narray_type rmnc_i = (p->rmnc())[s_cut];
		Rmnc_[i] = ifactory->interpolate_data( s_range, dblock_adapter(rmnc_i));
		narray_type zmns_i = (p->zmns())[s_cut];
		Zmns_[i] = ifactory->interpolate_data( s_range, dblock_adapter(zmns_i));
	}
}

morphism_vmec::~morphism_vmec() {
	if(Rmnc_) {
		for(size_t i = 0; i < xm_.size(); ++i)
			if(Rmnc_[i]) delete Rmnc_[i];
		delete Rmnc_;
	}
	if(Zmns_) {
		for(size_t i = 0; i < xm_.size(); ++i)
			if(Zmns_[i]) delete Zmns_[i];
		delete Zmns_;
	}
}

//! Maps `VMEC` coordinates @f$ q^\alpha @f$ into cartesian coordinates @f$ \textbf{x} @f$.
/*!
	In `VMEC` coordinates, @f$ \left( q^1, q^2, q^3 \right) = \left( s, \zeta, \theta \right) @f$
*/
IR3 morphism_vmec::operator()(const IR3& q) const {
	double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
	double R = 0.0, Z = 0.0;

	#pragma omp parallel for reduction(+: R, Z)
	for (size_t i = 0; i < xm_.size(); ++i) {  
		double m = xm_[i], n = xn_[i];
		double angle_mn = m*theta - n*zeta;
		double cosmn = std::cos( angle_mn );
		double sinmn = std::sin( angle_mn );
		// assuming for now that vmec equilibrium has stellarator symmetry.
		R += (*Rmnc_[i])(s) * cosmn; 
		Z += (*Zmns_[i])(s) * sinmn;
	}

	return {R * std::cos(zeta), R * std::sin(zeta), Z};
}

//! Inverse transform from cartesian coordinates @f$ \textbf{x} @f$ into `VMEC` coordinates @f$ q^\alpha @f$.
IR3 morphism_vmec::inverse(const IR3& X) const {
	typedef std::array<double, 2> IR2;
	double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
	double zeta = std::atan2(y, x);
	double r = std::sqrt(x*x + y*y);

	std::function<IR2(const IR2&)> zero_function =
		[&](const IR2& args) {
			auto [s, theta] = reflection_past_axis(args[0], args[1]);
			IR3 cyl = transform2cylindrical(IR3({s, zeta, theta}));
			double R = cyl[IR3::u], Z = cyl[IR3::w];
			return IR2({R-r, Z-z});
		};
	IR3 axis = transform2cylindrical({0, zeta, 0});
	double R0 = axis[IR3::u], Z0 = axis[IR3::w];
	IR2 guess = {0.5, std::atan2(z-Z0, r-R0)};
	IR2 roots = multiroot(1.0e-13, 100)(zero_function, guess);
	auto [s, theta] = reflection_past_axis(roots[0], roots[1]);
	return {s, zeta, theta};
}

//! Translates the curvilinear position `q` by the cartesian vector `delta`
/*!
	Implements the rule:
	@f$ q^\alpha = Q^\alpha \left( \textbf{X}\left(\right) + \Delta \right) @f$
*/
IR3 morphism_vmec::translation(const IR3 &q, const IR3 &delta) const {
	typedef std::array<double, 2> IR2;
	IR3 X = (*this)(q);
	double Xt = X[IR3::u] + delta[IR3::u];
	double Yt = X[IR3::v] + delta[IR3::v];
	double Zt = X[IR3::w] + delta[IR3::w];
	double Rt = std::sqrt(Xt*Xt + Yt*Yt);
	double zeta = std::atan2(Yt, Xt);
	std::function<IR2(const IR2&)> zero_function =
		[&](const IR2& args) {
			auto [s, theta] = reflection_past_axis(args[0], args[1]);
			IR3 cyl = transform2cylindrical(IR3({s, zeta, theta}));
			double R = cyl[IR3::u], Z = cyl[IR3::w];
			return IR2({R-Rt, Z-Zt});
		};
	IR2 guess = {q[IR3::u], q[IR3::w]};
	IR2 roots = multiroot(1.0e-12, 100)(zero_function, guess);
	auto [s, theta] = reflection_past_axis(roots[0], roots[1]);
	return {s, zeta, theta};
}

//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point @f$ q^\alpha @f$.
dIR3 morphism_vmec::del(const IR3& q) const {
	double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
	double R = 0.0, dR_ds = 0.0, dR_dtheta = 0.0, dR_dzeta = 0.0;
	double dZ_ds = 0.0, dZ_dtheta = 0.0, dZ_dzeta = 0.0;
	double sn_zeta = std::sin(zeta);
	double cn_zeta = std::cos(zeta);

	#pragma omp parallel for reduction(+: R, dR_ds, dR_dtheta, dR_dzeta, dZ_ds, dZ_dtheta, dZ_dzeta)
	for (size_t i = 0; i < xm_.size(); ++i) {  
		double m = xm_[i]; double n = xn_[i];
		double angle_mn = m*theta - n*zeta;
		double cosmn = std::cos(angle_mn);
		double sinmn = std::sin(angle_mn);
		double rmnc_i = (*Rmnc_[i])(s); 
		double zmns_i = (*Zmns_[i])(s);
		// assuming for now that vmec equilibrium has stellarator symmetry.
		R += rmnc_i * cosmn;
		dR_ds += (*Rmnc_[i]).derivative(s) * cosmn; 
		dR_dtheta += m * rmnc_i * sinmn; 
		dR_dzeta += n * rmnc_i * sinmn;
		dZ_ds += (*Zmns_[i]).derivative(s) * sinmn; 
		dZ_dtheta += m * zmns_i * cosmn; 
		dZ_dzeta += n * zmns_i * cosmn; 
	}
	dR_dtheta = -dR_dtheta;
	dZ_dzeta  = -dZ_dzeta;

	return {
		dR_ds*cn_zeta, dR_dzeta*cn_zeta - R*sn_zeta, dR_dtheta*cn_zeta,
		dR_ds*sn_zeta, dR_dzeta*sn_zeta + R*cn_zeta, dR_dtheta*sn_zeta,
		dZ_ds		 , dZ_dzeta					   , dZ_dtheta
	};
}

//! Returns the morphism's second derivatives, calculated in point `q`.
ddIR3 morphism_vmec::ddel(const IR3& q) const {
	double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
	double cn = std::cos(zeta);
	double sn = std::sin(zeta);

	double R = 0.0, dR_ds = 0.0, dR_dtheta = 0.0, dR_dzeta = 0.0;
	double d2R_ds2 = 0.0, d2R_dsdtheta = 0.0, d2R_dsdzeta = 0.0;
	double d2R_dzeta2 = 0.0, d2R_dzetadtheta = 0.0, d2R_dtheta2 = 0.0;

	double dZ_ds = 0.0, dZ_dtheta = 0.0, dZ_dzeta = 0.0;
	double d2Z_ds2 = 0.0, d2Z_dsdtheta = 0.0, d2Z_dsdzeta = 0.0;
	double d2Z_dzeta2 = 0.0, d2Z_dzetadtheta = 0.0, d2Z_dtheta2 = 0.0;

	#pragma omp parallel for reduction(+: R, dR_ds, dR_dtheta, dR_dzeta, d2R_ds2, d2R_dsdtheta, d2R_dsdzeta, d2R_dtheta2, d2R_dthetadzeta, d2R_dzeta2, Z, dZ_ds ,dZ_dtheta, dZ_dzeta, d2Z_ds2, d2Z_dsdtheta, d2Z_dsdzeta, d2Z_dtheta2, d2Z_dthetadzeta, d2Z_dzeta2)
	for (size_t i = 0; i < xm_.size(); ++i) {  
		double m = xm_[i]; double n = xn_[i];
		double cosmn = std::cos( m*theta - n*zeta );
		double sinmn = std::sin( m*theta - n*zeta );
		double rmnc_i = (*Rmnc_[i])(s); 
		double zmns_i = (*Zmns_[i])(s);
		double d_rmnc_i = (*Rmnc_[i]).derivative(s); 
		double d_zmns_i = (*Zmns_[i]).derivative(s); 
		double d2_rmnc_i = (*Rmnc_[i]).derivative2(s);
		double d2_zmns_i = (*Zmns_[i]).derivative2(s);
		// assuming for now that vmec equilibrium has stellarator symmetry.
		R += rmnc_i * cosmn;
		dR_ds += d_rmnc_i * cosmn; 
		dR_dzeta += n * rmnc_i * sinmn;
		dR_dtheta -= m * rmnc_i * sinmn;

		d2R_ds2 += d2_rmnc_i * cosmn; 
		d2R_dsdzeta += n * d_rmnc_i * sinmn;
		d2R_dsdtheta -= m * d_rmnc_i * sinmn;
		d2R_dzeta2 -= n * n * rmnc_i * cosmn;
		d2R_dzetadtheta += m * n * rmnc_i * cosmn;
		d2R_dtheta2 -= m * m * rmnc_i * cosmn;

		dZ_ds += d_zmns_i * sinmn; 
		dZ_dzeta -= n * zmns_i * cosmn; 
		dZ_dtheta += m * zmns_i * cosmn; 

		d2Z_ds2 += d2_zmns_i * sinmn;
		d2Z_dsdzeta -= n * d_zmns_i * cosmn;
		d2Z_dsdtheta += m * d_zmns_i * cosmn;
		d2Z_dzeta2 -= n * n * zmns_i * sinmn;
		d2Z_dzetadtheta += m * n * zmns_i * sinmn;
		d2Z_dtheta2 -= m * m * zmns_i * sinmn;
	}

	return {
		d2R_ds2*cn, d2R_dsdzeta*cn-dR_ds*sn, d2R_dsdtheta*cn, (d2R_dzeta2-R)*cn-2*dR_dzeta*sn, d2R_dzetadtheta*cn-dR_dtheta*sn, d2R_dtheta2*cn,
		d2R_ds2*sn, d2R_dsdzeta*sn+dR_ds*cn, d2R_dsdtheta*sn, (d2R_dzeta2-R)*sn+2*dR_dzeta*cn, d2R_dzetadtheta*sn+dR_dtheta*cn, d2R_dtheta2*sn,
		d2Z_ds2   , d2Z_dsdzeta            , d2Z_dsdtheta   ,  d2Z_dzeta2                    , d2Z_dzetadtheta                , d2Z_dtheta2
	};
}

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
double morphism_vmec::jacobian(const IR3&q) const {
	double s = q[IR3::u], zeta = q[IR3::v], theta = q[IR3::w];
	double R = 0.0, dR_ds = 0.0, dR_dtheta = 0.0;
	double dZ_ds = 0.0, dZ_dtheta = 0.0;
	double sn_zeta = std::sin(zeta);
	double cn_zeta = std::cos(zeta);

	#pragma omp parallel for reduction(+: R, Z, dR_ds, dR_dtheta, dR_dzeta, dZ_ds, dZ_dtheta, dZ_dzeta)
	for (size_t i = 0; i < xm_.size(); ++i) {  
		double m = xm_[i]; double n = xn_[i];
		double angle_mn = m*theta - n*zeta;
		double cosmn = std::cos(angle_mn);
		double sinmn = std::sin(angle_mn);
		double rmnc_i = (*Rmnc_[i])(s); 
		double zmns_i = (*Zmns_[i])(s);
		// assuming for now that vmec equilibrium has stellarator symmetry.
		R += rmnc_i * cosmn;
		dR_ds += (*Rmnc_[i]).derivative(s) * cosmn; 
		dR_dtheta += m * rmnc_i * sinmn;
		dZ_ds += (*Zmns_[i]).derivative(s) * sinmn; 
		dZ_dtheta += m * zmns_i * cosmn;
	}
	dR_dtheta = -dR_dtheta;

	return R * (dR_ds * dZ_dtheta - dR_dtheta * dZ_ds);
}

// //@todo move this to jacobian and think about testing this by calling the parent
// double morphism_vmec::jacobian_vmec(const IR3& position) const {
// 	double s = position[IR3::u];
// 	double zeta = position[IR3::v];
// 	double theta = position[IR3::w];
// 	double J = 0.0;
// 	#pragma omp parallel for reduction(+: J)
// 	for (size_t i = 0; i < xm_nyq_.size(); i++) {  
// 		J += (*gmnc_[i])(s) * std::cos( xm_nyq_[i]*theta - xn_nyq_[i]*zeta );
// 	};
// 	// left-handed VMEC coordinate system is re-oriented 
// 	// to u = Phi/Phi_bnd, v = zeta, w = theta for J>0
// 	// should we check/assume that signgs is always negative?
// 	return -J;
// }

IR3 morphism_vmec::transform2cylindrical(const IR3& position) const {
	double u = position[gyronimo::IR3::u];
	double v = position[gyronimo::IR3::v];
	double w = position[gyronimo::IR3::w];
	double R = 0.0, Z = 0.0;

	#pragma omp parallel for reduction(+: R, Z)
	for (size_t i = 0; i < xm_.size(); ++i) {
		double m = xm_[i]; double n = xn_[i];
		R+= (*Rmnc_[i])(u) * std::cos( m*w - n*v ); 
		Z+= (*Zmns_[i])(u) * std::sin( m*w - n*v );
	}
	return  {R, v, Z};
}

// double morphism_vmec::reduce_2pi(double angle) const {
// 	double l = 2*std::numbers::pi;
// 	angle -= l*std::floor(angle/l);
// 	return angle;
// }

std::tuple<double, double> morphism_vmec::reflection_past_axis(
		double s, double theta) const {
	if(s < 0)
		// return {-s, reduce_2pi(theta + std::numbers::pi)};
		return {-s, theta + std::numbers::pi};
	else
		return {s, theta};
}

} // end namespace gyronimo
