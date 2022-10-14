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

// @lorentz.cc, this file is part of ::gyronimo::

#include <gyronimo/dynamics/lorentz.hh>
#include <cmath>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>

namespace gyronimo {

//! Class constructor.
/*!
	Takes a reference length `Lref`, a reference velocity `Vref`, 
	the charge-over-mass ratio `qom`, an electric field `Efield`
	and a magnetic field `Bfield`. The fields `Efield` and `Bfield`
	cannot be `nullptr`. For null fields, create an IR3field that
	returns zero in the contravariant components.
*/
lorentz::lorentz(const double &Lref, const double &Vref, 
		const double &qom, const IR3field *Efield, const IR3field *Bfield)
		: Lref_(Lref), Vref_(Vref), Tref_(Lref/Vref), 
		Oref_(Bfield ? qom * codata::e / 
			codata::m_proton * Bfield->m_factor() * Tref_ : 1.0), 
		electric_field_(Efield), magnetic_field_(Bfield),
		iEfield_time_factor_(Efield ? Tref_ / Efield->t_factor() : 1.0),
		iBfield_time_factor_(Bfield ? Tref_ / Bfield->t_factor() : 1.0),
		metric_(Bfield ? Bfield->metric() : nullptr) {

	// test if fields exist
	// if(!Efield) error(__func__, __FILE__, __LINE__, 
	//             " Efield cannot be nullptr.", 1);
	if(!Bfield) error(__func__, __FILE__, __LINE__, 
	            " Bfield cannot be nullptr.", 1);

	if(Efield) {
		// test if fields are in the same coordinates
		const metric_covariant *Emetric = Efield->metric();
		const metric_covariant *Bmetric = Bfield->metric();
		if(Emetric != Bmetric)
			error(__func__, __FILE__, __LINE__, 
				" Efield and Bfield must be in the same coordinates.", 1);

		// test if normalization factors are consistent
		double E_m_factor = Efield->m_factor();
		double B_m_factor = Bfield->m_factor();
		if(std::abs(1.0 - E_m_factor/(Vref_*B_m_factor)) > 1.0e-14)
			error(__func__, __FILE__, __LINE__, 
				" inconsistent Efield and Bfield.", 1);
	}
}

//! Calculates the dynamic system `dx` for a particle in state `state`, at time `t`.
lorentz::state lorentz::operator()(const state &s, const double &time) const {

	// extract position and velocity
	IR3 q = {s[0], s[1], s[2]};
	IR3 v = {s[3], s[4], s[5]};
	
	IR3 dq = Lref_ * v;
	IR3 dv = Lref_ * metric_->inertial_force(q, v);

	if(electric_field_) {
		double Etime = time * iEfield_time_factor_;
		dv += Oref_ * electric_field_->contravariant(q, Etime);
	}
	if(magnetic_field_) {
		double Btime = time * iBfield_time_factor_;
		IR3 B = magnetic_field_->contravariant(q, Btime);
		double jacobian = metric_->jacobian(q);
		IR3 vxB_covariant = cross_product<covariant>(v, B, jacobian);
		IR3 vxB = metric_->to_contravariant(vxB_covariant, q);
		dv += Oref_ * vxB;
	}

	return {dq[IR3::u], dq[IR3::v], dq[IR3::w],
	        dv[IR3::u], dv[IR3::v], dv[IR3::w]};
}

//! Generates a `lorentz::state` from curvilinear position `pos` and velocity `vel`.
lorentz::state lorentz::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

//! Extracts the position from a `lorentz::state`.
IR3 lorentz::get_position(const state &s) const {
	return {s[0], s[1], s[2]};
}

//! Extracts the velocity from a `lorentz::state`, normalized to `Vref`.
IR3 lorentz::get_velocity(const state &s) const {
	return {s[3], s[4], s[5]};
}

//! Returns the kinetic energy of the state, normalized to `Uref`.
double lorentz::energy_kinetic(const state &s) const {
	IR3 pos   = {s[0], s[1], s[2]};
	IR3 v_con = {s[3], s[4], s[5]};
	IR3 v_cov = metric_->to_covariant(v_con, pos);
	return inner_product(v_cov, v_con);
}

//! Returns the parallel energy of the state, normalized to `Uref`.
double lorentz::energy_parallel(const state &s, const double &time) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 v_con = {s[3], s[4], s[5]};
	double Btime = time * iBfield_time_factor_;
	IR3 b_cov = magnetic_field_->covariant_versor(q, Btime);
	double vpar = inner_product(v_con, b_cov);
	return vpar*vpar;
}

//! Returns the perpendicular energy of the state, normalized to `Uref`.
double lorentz::energy_perpendicular(const state &s, const double &time) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 v_con = {s[3], s[4], s[5]};
	double Btime = time * iBfield_time_factor_;
	IR3 b_con = magnetic_field_->contravariant_versor(q, Btime);
	double jacobian = metric_->jacobian(q);
	IR3 vperp_cov = cross_product<covariant>(v_con, b_con, jacobian);
	IR3 vperp_con = metric_->to_contravariant(vperp_cov, q);
	return inner_product(vperp_cov, vperp_con);
}

} // end namespace gyronimo