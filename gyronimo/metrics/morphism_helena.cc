// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues and Manuel Assunção.

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

// @morphism_helena.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_helena.hh>
#include <numbers>
#include <gyronimo/core/multiroot.hh>

namespace gyronimo{

morphism_helena::morphism_helena(
		const parser_helena *p, const interpolator2d_factory *ifactory)
		: parser_(p), R_(nullptr), z_(nullptr) {
	double Rgeo = p->rgeo();
	double a = p->eps()*Rgeo;
	dblock_adapter s_range(p->s()), chi_range(p->chi());
	R_ = ifactory->interpolate_data(
		s_range, chi_range, dblock_adapter(
			parser_helena::narray_type(a*p->x() + Rgeo)));
	z_ = ifactory->interpolate_data(
		s_range, chi_range, dblock_adapter(
			parser_helena::narray_type(a*p->y() + 0.0)));
}

morphism_helena::~morphism_helena() {
	if(R_) delete R_;
	if(z_) delete z_;
}

//! Maps `HELENA` coordinates @f$ q^\alpha @f$ into cartesian coordinates @f$ \textbf{x} @f$.
/*!
	In `HELENA` coordinates, @f$ \left( q^1, q^2, q^3 \right) = \left( s, \chi, \phi \right) @f$,
	where @f$ s @f$ is the flux surface label, @f$ \chi @f$ is the poloidal angle 
	and @f$ \phi @f$ is the toroidal angle.
	Implements the coordinate transformation:
	@f$ \left( x, y, z \right) = \textbf{R} \left( s, \chi, \phi \right) @f$
*/
IR3 morphism_helena::operator()(const IR3& q) const {
	double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]), phi = q[IR3::w];
	double R = (*R_)(s, chi);
	return {R*std::cos(phi), -R*std::sin(phi), (*z_)(s, chi)};
}

//! Inverse transform from cartesian coordinates @f$ \textbf{x} @f$ into `HELENA` coordinates @f$ q^\alpha @f$.
/*!
	Implements the inverse transformation:
	@f$ \left( s, \chi, \phi \right) = \textbf{R}^{-1} \left( x, y, z \right) @f$
*/
IR3 morphism_helena::inverse(const IR3& X) const {
	typedef std::array<double, 2> IR2;
	double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
	double R = std::sqrt(x*x + y*y);
	std::function<IR2(const IR2&)> zero_function =
		[&](const IR2& args) {
			auto [s, chi] = reflection_past_axis(args[0], args[1]);
			return IR2({(*R_)(s, chi) - R, (*z_)(s, chi) - z});
		};
	IR2 guess = {0.5, std::atan2(z, R - parser_->rmag())};
	IR2 roots = multiroot(1.0e-15, 100)(zero_function, guess);
	auto [s, chi] = reflection_past_axis(roots[0], roots[1]);
	return {s, chi, std::atan2(-y, x)};
}

//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point @f$ q^\alpha @f$.
dIR3 morphism_helena::del(const IR3& q) const {
	double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]), phi = q[IR3::w];
	double R = (*R_)(s, chi),
		Ru = (*R_).partial_u(s, chi), Rv = (*R_).partial_v(s, chi);
	double cos = std::cos(phi), sin = std::sin(phi);
	return {Ru*cos, Rv*cos, -R*sin, -Ru*sin, -Rv*sin, -R*cos,
		(*z_).partial_u(s, chi), (*z_).partial_v(s, chi), 0.0};
}

//! Returns the morphism's second derivatives, calculated in point @f$ q^\alpha @f$.
ddIR3 morphism_helena::ddel(const IR3 &q) const {
	double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]), phi = q[IR3::w];
	double R = (*R_)(s, chi), Ru = R_->partial_u(s, chi), Rv = R_->partial_v(s, chi);
	double Ruu = R_->partial2_uu(s, chi), Ruv = R_->partial2_uv(s, chi), Rvv = R_->partial2_vv(s, chi);
	double Zu = z_->partial_u(s, chi), Zv = z_->partial_v(s, chi);
	double Zuu = z_->partial2_uu(s, chi), Zuv = z_->partial2_uv(s, chi), Zvv = z_->partial2_vv(s, chi);
	double sn = std::sin(phi), cn = std::cos(phi);
	return {
	//	   iuu      iuv     iuw      ivv     ivw    iww
		 Ruu*cn,  Ruv*cn, -Ru*sn,  Rvv*cn, -Rv*sn, -R*cn,
		-Ruu*sn, -Ruv*sn, -Ru*cn, -Rvv*sn, -Rv*cn,  R*sn,
		    Zuu,  Zuv   ,    0  ,  Zvv   ,    0  ,   0
	};
}

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
double morphism_helena::jacobian(const IR3 &q) const {
	double s = q[IR3::u], chi = parser_->reduce_chi(q[IR3::v]);
	double R = (*R_)(s, chi), Ru = R_->partial_u(s, chi), Rv = R_->partial_v(s, chi);
	double Zu = z_->partial_u(s, chi), Zv = z_->partial_v(s, chi);
	return R * (Ru * Zv - Rv * Zu);
}

std::tuple<double, double> morphism_helena::reflection_past_axis(
		double s, double chi) const {
	if(s < 0)
		return {-s, parser_->reduce_chi(chi + std::numbers::pi)};
	else
		return {s, parser_->reduce_chi(chi)};
}

} // end namespace gyronimo
