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

// @morphism_polar_torus.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/morphism_polar_torus.hh>
#include <cmath>

namespace gyronimo {

morphism_polar_torus::morphism_polar_torus(const double minor_radius, double major_radius)
    : minor_radius_(minor_radius), major_radius_(major_radius),
      iaspect_ratio_(minor_radius_/major_radius_),
	  volume_factor_(minor_radius*minor_radius*major_radius),
	  iminor_radius_(1.0/minor_radius) {
}

//! Maps polar torus coordinates @f$ q^\alpha @f$ into cartesian coordinates @f$ \textbf{x} @f$.
/*!
	In polar torus coordinates, @f$ \left( q^1, q^2, q^3 \right) = \left( r, \theta, \phi \right) @f$.
	Implements the coordinate transformation:
	@f{gather*}{
		\begin{aligned}
			x &=  \left( R_{0} + a \, r \, \cos \theta \right) \, \cos \phi \\
			y &= -\left( R_{0} + a \, r \, \cos \theta \right) \, \sin \phi \\
			z &= a \, r \, \sin \theta
		\end{aligned}
    @f}
*/
IR3 morphism_polar_torus::operator()(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::v]), sn_theta = std::sin(q[IR3::v]);
	double cn_phi = std::cos(q[IR3::w]), sn_phi = std::sin(q[IR3::w]);
	double R = major_radius_*(1.0 + iaspect_ratio_*r*cn_theta);
	return { R*cn_phi, 
			-R*sn_phi, 
			minor_radius_*r*sn_theta};
}

//! Inverse transform from cartesian coordinates @f$ \textbf{x} @f$ into spherical coordinates @f$ q^\alpha @f$.
/*!
	Implements the inverse transformation:
	@f{gather*}{
		\begin{aligned}
			r &= \frac{1}{a} \, \sqrt{z^2 + (\sqrt{x^2 + y^2} - R_{0})^2} \\
			\theta &= \arctan \left( \frac{z}{\sqrt{x^2 + y^2} - R_{0}} \right) \\
			\phi &= -\arctan \left( \frac{y}{x} \right)
		\end{aligned}
    @f}
*/
IR3 morphism_polar_torus::inverse(const IR3 &x) const {
	double z = x[IR3::w];
	double R = std::sqrt(x[IR3::u]*x[IR3::u] + x[IR3::v]*x[IR3::v]);
	double dR = R - major_radius_;
	return {iminor_radius_*std::sqrt(z*z + dR*dR), 
			std::atan2(z, dR), 
			std::atan2(-x[IR3::v], x[IR3::u])};
}

//! Returns the morphism's first derivatives, correspondent to the covariant basis vectors in point @f$ q^\alpha @f$.
/*!
	Implements the coordinate transformation's first derivatives:
	@f{gather*}{
		\begin{aligned}
			\textbf{e}_r &= \left( a \, \cos \theta \, \cos \phi, -a \, \cos \theta \, \sin \phi, a \, \sin \theta \right) \\
			\textbf{e}_\theta &= \left( -a \, r \, \sin \theta \, \cos \phi, a \, r \, \sin \theta \, \sin \phi, a \, r \, \cos \theta \right) \\
			\textbf{e}_\phi &= \left( -\left( R_0 + a \, r \, \cos \theta \right) \, \sin \phi, \left( R_0 + a \, r \, \cos \theta \right) \, \cos \phi, 0 \right)
		\end{aligned}
    @f}
*/
dIR3 morphism_polar_torus::del(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::v]);
	double sn_theta = std::sin(q[IR3::v]);
	double cn_phi = std::cos(q[IR3::w]);
	double sn_phi = std::sin(q[IR3::w]);
	double R = major_radius_*(1.0 + iaspect_ratio_*r*cn_theta);
	double a_cn = minor_radius_*cn_theta;
	double ar_cn = r * a_cn;
	double a_sn = minor_radius_*sn_theta;
	double ar_sn = r * a_sn;
	return {
		 a_cn * cn_phi, -ar_sn * cn_phi, -R * sn_phi,
		-a_cn * sn_phi,  ar_sn * sn_phi, -R * cn_phi,
		 a_sn         ,  ar_cn         ,  0
	};
}

//! Returns the morphism's second derivatives, calculated in point @f$ q^\alpha @f$.
ddIR3 morphism_polar_torus::ddel(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::v]);
	double sn_theta = std::sin(q[IR3::v]);
	double cn_phi = std::cos(q[IR3::w]);
	double sn_phi = std::sin(q[IR3::w]);
	double R = major_radius_*(1.0 + iaspect_ratio_*r*cn_theta);
	double a_cn = minor_radius_*cn_theta;
	double ar_cn = r * a_cn;
	double a_sn = minor_radius_*sn_theta;
	double ar_sn = r * a_sn;
	return {
	// iuu      iuv             iuw              ivv             ivw          iww
		0, -a_sn * cn_phi, -a_cn * sn_phi, -ar_cn * cn_phi, ar_sn * sn_phi, -R * cn_phi,
		0,  a_sn * sn_phi, -a_cn * cn_phi,  ar_cn * sn_phi, ar_sn * cn_phi,  R * sn_phi,
		0,  a_cn         ,       0       , -ar_sn         ,       0       ,  0
	};
}

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
/*!
	Implements the Jacobian in spherical coordinates: 
	@f$ J = a^2 \, r \, \left( R_0 + a \, r \, \cos \theta \right) @f$
*/
double morphism_polar_torus::jacobian(const IR3 &q) const {
	double r = q[IR3::u];
	double cn = std::cos(q[IR3::v]);
	double Rfactor = (1.0+iaspect_ratio_*r*cn);
	return volume_factor_*r*Rfactor;
}

//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point @f$ q^\alpha @f$.
/*!
	Implements the inverse transformation's first derivatives:
	@f{gather*}{
		\begin{aligned}
			\textbf{e}^r &= \left( \frac{\cos \theta \, \cos \phi}{a}, -\frac{\cos \theta \, \sin \phi}{a}, \frac{\sin \theta}{a} \right) \\
			\textbf{e}^\theta &= \left( -\frac{\sin \theta \, \cos \phi}{a \, r}, \frac{\sin \theta \, \sin \phi}{a \, r}, \frac{\cos \theta}{a \, r} \right) \\
			\textbf{e}^\phi &= \left( - \frac{\sin \phi}{R_0 + a \, r \, \cos \theta}, \frac{\cos \phi}{R_0 + a \, r \, \cos \theta}, 0 \right)
		\end{aligned}
    @f}
*/
dIR3 morphism_polar_torus::del_inverse(const IR3 &q) const {
	double r = q[IR3::u];
	double cn_theta = std::cos(q[IR3::v]);
	double sn_theta = std::sin(q[IR3::v]);
	double cn_phi = std::cos(q[IR3::w]);
	double sn_phi = std::sin(q[IR3::w]);
	double R = major_radius_*(1.0 + iaspect_ratio_*r*cn_theta);
	double iR = 1.0 / R;
	double ia_cn = iminor_radius_*cn_theta;
	double iar_cn = r * ia_cn;
	double ia_sn = iminor_radius_*sn_theta;
	double iar_sn = r * ia_sn;
	return {
		  ia_cn * cn_phi, -ia_cn * sn_phi,  ia_sn,
		-iar_sn * cn_phi, iar_sn * sn_phi, iar_cn,
		    -iR * sn_phi,    -iR * cn_phi, 0
	};
}

} // end namespace gyronimo