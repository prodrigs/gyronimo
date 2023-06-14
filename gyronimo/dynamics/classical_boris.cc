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

// @classical_boris.cc, this file is part of ::gyronimo::

#include <gyronimo/dynamics/classical_boris.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/core/error.hh>
#include <gyronimo/dynamics/boris_push.hh>
#include <cmath>

#include <gyronimo/dynamics/lorentz.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>
#include <boost/numeric/odeint.hpp>

namespace gyronimo {

/*!
	Class constructor.
	Takes a reference length `Lref`, a reference velocity `Vref`, 
	the charge-over-mass ratio `qom`, an electric field `Efield`
	and a magnetic field `Bfield`. The field `Bfield`
	cannot be `nullptr`.
*/
classical_boris::classical_boris(const double &Lref, const double &Vref, 
		const double &qom, const IR3field *Efield, const IR3field *Bfield)
		: Lref_(Lref), Vref_(Vref), Tref_(Lref/Vref), qom_(qom),
		Oref_(Bfield ? qom * codata::e / 
			codata::m_proton * Bfield->m_factor() * Tref_ : 1.0), 
		electric_field_(Efield), magnetic_field_(Bfield),
		iEfield_time_factor_(Efield ? Tref_ / Efield->t_factor() : 1.0),
		iBfield_time_factor_(Bfield ? Tref_ / Bfield->t_factor() : 1.0),
		metric_(Bfield ? dynamic_cast<const metric_connected*>(Bfield->metric()) : nullptr),
		my_morphism_(metric_ ? metric_->my_morphism() : nullptr) {

	// test if fields exist
	// if(!Efield) error(__func__, __FILE__, __LINE__, 
	// 	    " Efield cannot be nullptr.", 1);
	if(!Bfield) error(__func__, __FILE__, __LINE__, 
		    " Bfield cannot be nullptr.", 1);

	// test if metrics are the same
	if(Efield) {
		const metric_covariant *Emetric = Efield->metric();
		const metric_covariant *Bmetric = Bfield->metric();
		if(Emetric != Bmetric) error(__func__, __FILE__, __LINE__, 
			" Efield and Bfield must be in the same coordinate representation.", 1);
	}

	// ensure that the field metrics derive from metric_connected
	if(!metric_) error(__func__, __FILE__, __LINE__, 
		" field metrics must derive from 'gyronimo::metric_connected'.", 1);

	// test if normalization factors are consistent
	if(Efield) {
		double E_m_factor = Efield->m_factor();
		double B_m_factor = Bfield->m_factor();
		if(std::abs(1.0 - E_m_factor/(Vref_*B_m_factor)) > 1.0e-14)
			error(__func__, __FILE__, __LINE__, 
				" inconsistent Efield and Bfield normalization.", 1);
	}
}

//! Performs a time step `dt` update to the state `in` and returns the result.
classical_boris::state classical_boris::do_step(
		const state &in, const double &time, const double &dt) const {
	
	// extract position and velocity from state
	IR3 q_old = {in[0], in[1], in[2]};
	IR3 v_old = {in[3], in[4], in[5]};
	dIR3 e_old = my_morphism_->del(q_old);

	// calculate fields
	IR3 Efield = {0.0, 0.0, 0.0};
	if(electric_field_) {
		double Etime = time * iEfield_time_factor_;
		IR3 E_contravariant = electric_field_->contravariant(q_old, Etime);
		Efield = contraction<second>(e_old, E_contravariant);
	}
	double Btime = time * iBfield_time_factor_;
	IR3 Bfield_contravariant = magnetic_field_->contravariant(q_old, Btime);
	IR3 Bfield = contraction<second>(e_old, Bfield_contravariant);
	double Bmag = std::sqrt(inner_product(Bfield, Bfield));
    IR3 Bversor = Bfield / Bmag;

	// perform cartesian_boris step
	IR3 v_new = electric_field_ ? 
		boris_push(v_old, Oref_, Efield, Bmag, Bversor, dt) :
		boris_push(v_old, Oref_, Bmag, Bversor, dt);
	// IR3 x_new = x_old + (Lref_*dt) * v_new;

	// IR3 q_new = my_morphism_->inverse(x_new);
	// IR3 u_new = my_morphism_->to_contravariant(v_new, q_new);

	IR3 q_new = my_morphism_->translation(q_old, (Lref_*dt) * v_new);

	return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w],
			v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};
}

//! Returns the `classical_boris::state` state from a point in SI phase-space.
classical_boris::state classical_boris::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

//! Returns the vector position of the state in SI units.
IR3 classical_boris::get_position(const state &s) const {
	return {s[0], s[1], s[2]};
}

//! Returns the vector velocity of the state in SI units.
IR3 classical_boris::get_velocity(const state &s) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 v = {s[3], s[4], s[5]};
	return my_morphism_->to_contravariant(v, q);
}

//! Returns the kinetic energy of the state, normalized to `Uref`.
double classical_boris::energy_kinetic(const state &s) const {
	return (s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
}

//! Returns the parallel energy of the state, normalized to `Uref`.
double classical_boris::energy_parallel(const state &s, double &time) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 v = {s[3], s[4], s[5]};
	double Btime = time * iBfield_time_factor_;
	IR3 Bversor_con = magnetic_field_->contravariant_versor(q, Btime);
	IR3 Bversor = my_morphism_->from_contravariant(Bversor_con, q);
	double vpar = inner_product(v, Bversor);
	return vpar*vpar;
}

//! Returns the perpendicular energy of the state, normalized to `Uref`.
double classical_boris::energy_perpendicular(const state &s, double &time) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 v = {s[3], s[4], s[5]};
	double Btime = time * iBfield_time_factor_;
	IR3 Bversor_con = magnetic_field_->contravariant_versor(q, Btime);
	IR3 Bversor = my_morphism_->from_contravariant(Bversor_con, q);
	IR3 vperp = cross_product(v, Bversor);
	return inner_product(vperp, vperp);
}

classical_boris::state classical_boris::generate_initial_state(
		const IR3 &pos, const IR3 &vel, 
		const double &time, const double &dt) const {

	// state s0 = {
	// 	pos[IR3::u], pos[IR3::v], pos[IR3::w],
	// 	vel[IR3::u], vel[IR3::v], vel[IR3::w]
	// };
	// state s1 = do_step(s0, time, -0.5*dt);

	// return {s0[0], s0[1], s0[2], s1[3], s1[4], s1[5]};
    
    lorentz lo(Lref_, Vref_, qom_, electric_field_, magnetic_field_);
    odeint_adapter<gyronimo::lorentz> sys(&lo);
    boost::numeric::odeint::runge_kutta4<gyronimo::lorentz::state> rk4;
    IR3 vel_con = my_morphism_->to_contravariant(vel, pos);
	lorentz::state inout = {
		pos[IR3::u], pos[IR3::v], pos[IR3::w],
		vel_con[IR3::u], vel_con[IR3::v], vel_con[IR3::w]
	};
	rk4.do_step(sys, inout, time, -0.5*dt);
    IR3 qmh = {inout[0], inout[1], inout[2]};
    IR3 vmh_con = {inout[3], inout[4], inout[5]};
    IR3 vmh = my_morphism_->from_contravariant(vmh_con, qmh);

	return {pos[IR3::u], pos[IR3::v], pos[IR3::w], 
            vmh[IR3::u], vmh[IR3::v], vmh[IR3::w]};
}

} // end namespace gyronimo