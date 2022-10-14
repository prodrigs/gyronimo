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

// @cartesian_boris.cc, this file is part of ::gyronimo::

#include <gyronimo/dynamics/cartesian_boris.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/core/error.hh>
#include <gyronimo/metrics/metric_cartesian.hh>
#include <gyronimo/dynamics/boris_push.hh>

namespace gyronimo {

/*!
	Class constructor.
	Takes a reference length `Lref`, a reference velocity `Vref`, 
	the charge-over-mass ratio `qom`, an electric field `Efield`
	and a magnetic field `Bfield`. The field `Bfield`
	cannot be `nullptr`.
*/
cartesian_boris::cartesian_boris(const double &Lref, const double &Vref, 
		const double &qom, const IR3field *Efield, const IR3field *Bfield)
		: Lref_(Lref), Vref_(Vref), Tref_(Lref/Vref), qom_(qom),
		Oref_(Bfield ? qom * codata::e / 
			codata::m_proton * Bfield->m_factor() * Tref_ : 1.0), 
		electric_field_(Efield), magnetic_field_(Bfield),
		iEfield_time_factor_(Efield ? Tref_ / Efield->t_factor() : 1.0),
		iBfield_time_factor_(Bfield ? Tref_ / Bfield->t_factor() : 1.0),
		metric_(Bfield ? dynamic_cast<const metric_cartesian*>(Bfield->metric()) : nullptr),
		morph_(metric_ ? metric_->morph() : nullptr) {

	// test if fields exist
	// if(!Efield) error(__func__, __FILE__, __LINE__, 
	// 							" Efield cannot be nullptr.", 1);
	if(!Bfield) error(__func__, __FILE__, __LINE__, 
								" Bfield cannot be nullptr.", 1);

	// test if fields are cartesian
	if(!metric_)
		error(__func__, __FILE__, __LINE__, 
			" Bfield must be in cartesian coordinates.", 1);
	if(Efield) {
		const metric_covariant *Emetric = Efield->metric();
		if(Emetric != metric_)
			error(__func__, __FILE__, __LINE__, 
				" Efield must be in cartesian coordinates.", 1);
	}

	// test if normalization factors are consistent
	if(Efield) {
		double E_m_factor = Efield->m_factor();
		double B_m_factor = Bfield->m_factor();
		if(std::abs(1.0 - E_m_factor/(Vref_*B_m_factor)) > 1.0e-14)
			error(__func__, __FILE__, __LINE__, 
				" inconsistent Efield and Bfield.", 1);
	}
}

//! Performs a time step `dt` update to the state `in` and returns the result.
cartesian_boris::state cartesian_boris::do_step(
		const cartesian_boris::state &in, const double &time, const double &dt) const {
	
	// extract position and velocity from state
	IR3 x_old = {in[0], in[1], in[2]};
	IR3 v_old = {in[3], in[4], in[5]};

	// calculate fields
	IR3 Efield = {0.0, 0.0, 0.0};
	if(electric_field_) {
		double Etime = time * iEfield_time_factor_;
		Efield = electric_field_->contravariant(x_old, Etime);
	}
	double Btime = time * iBfield_time_factor_;
	double Bmag = magnetic_field_->magnitude(x_old, Btime);
	IR3 Bversor = magnetic_field_->contravariant_versor(x_old, Btime);

	// perform cartesian_boris step
	IR3 v_new = electric_field_ ? 
		boris_push(v_old, Oref_, Efield, Bmag, Bversor, dt) :
		boris_push(v_old, Oref_, Bmag, Bversor, dt);
	IR3 x_new = x_old + (Lref_*dt) * v_new;

	return {x_new[IR3::u], x_new[IR3::v], x_new[IR3::w],
			v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};
}

//! Returns the `cartesian_boris::state` state from a point in SI phase-space.
cartesian_boris::state cartesian_boris::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

//! Returns the vector position of the state in SI units.
IR3 cartesian_boris::get_position(const state &s) const {
	return {s[0], s[1], s[2]};
}

//! Returns the vector velocity of the state in SI units.
IR3 cartesian_boris::get_velocity(const state &s) const {
	return {s[3], s[4], s[5]};
}

//! Returns the kinetic energy of the state, normalized to `Uref`.
double cartesian_boris::energy_kinetic(const state &s) const {
	return s[3]*s[3]+s[4]*s[4]+s[5]*s[5];
}

//! Returns the parallel energy of the state, normalized to `Uref`.
double cartesian_boris::energy_parallel(const state &s, double &time) const {
	IR3 x = {s[0], s[1], s[2]};
	IR3 v = {s[3], s[4], s[5]};
	double Btime = time * iBfield_time_factor_;
	IR3 b = magnetic_field_->contravariant_versor(x, Btime);
	double vpar = inner_product(v, b);
	return vpar*vpar;
}

//! Returns the perpendicular energy of the state, normalized to `Uref`.
double cartesian_boris::energy_perpendicular(const state &s, double &time) const {
	IR3 x = {s[0], s[1], s[2]};
	IR3 v = {s[3], s[4], s[5]};
	double Btime = time * iBfield_time_factor_;
	IR3 b = magnetic_field_->contravariant_versor(x, Btime);
	IR3 vperp = cross_product(v, b);
	return inner_product(vperp, vperp);
}

//! Creates the first `cartesian_boris::state` from a point in cartesian phase-space.
cartesian_boris::state cartesian_boris::generate_initial_state(
		const IR3 &cartesian_position, const IR3 &cartesian_velocity, 
		const double &time, const double &dt) const {

	state s0 = {
		cartesian_position[IR3::u],
		cartesian_position[IR3::v],
		cartesian_position[IR3::w],
		cartesian_velocity[IR3::u],
		cartesian_velocity[IR3::v],
		cartesian_velocity[IR3::w]
	};
	state s1 = do_step(s0, time, -0.5*dt);

	return {s0[0], s0[1], s0[2], s1[3], s1[4], s1[5]};
}

} // end namespace gyronimo