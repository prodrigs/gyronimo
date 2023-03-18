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

// @curvilinear_boris.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_CURVILINEAR_BORIS
#define GYRONIMO_CURVILINEAR_BORIS

#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/metrics/metric_connected.hh>

namespace gyronimo {

//! Gyronimo implementation for curvilinear boris stepper class.
/*!
	This class implements the curvilinear boris pusher 
	[Journal of Plasma Physics, Volume 61, Issue 3, April 1999, pp. 367 - 389].
	This class mixes the typical Cartesian boris push with a Runge-Kutta 2 advance
	in curvilinear coordinates. To achieve this, field representation is only 
	constrained to metrics of the type `metric_connected` to be able to convert 
	coordinates and vectors using the connected `morphism`.
	To use the classical boris pusher in curvilinear coordinates, look at
	`classical_boris` class.
*/
class curvilinear_boris {
public:
	using state = std::array<double,6>;

	curvilinear_boris(const double &Lref, const double &Vref, 
		const double &qom, const IR3field *Efield, const IR3field *Bfield);
	~curvilinear_boris() {};

	state do_step(const state &in, const double &time, const double &dt) const;

	state generate_state(const IR3 &pos, const IR3 &vel) const;
	IR3 get_position(const state &s) const;
	IR3 get_velocity(const state &s) const;

	double energy_kinetic(const state &s) const;
	double energy_parallel(const state &s, double &time) const;
	double energy_perpendicular(const state &s, double &time) const;

	state generate_initial_state(const IR3 &pos, const IR3 &vel, 
		const double &t, const double &dt) const;

	const IR3field* electric_field() const {return electric_field_;};
	const IR3field* magnetic_field() const {return magnetic_field_;};
	const morphism* my_morphism() const {return my_morphism_;};

	//! Returns the reference length scale `Lref`.
	double Lref() const {return Lref_;};
	//! Returns the reference time scale `Tref`.
	double Tref() const {return Tref_;};
	//! Returns the reference velocity scale `Vref`.
	double Vref() const {return Vref_;};
	//! Returns the reference frequency scale `Oref`.
	double Oref() const {return Oref_;};
	//! Returns the charge-over-mass ratio.
	double qom() const {return qom_;};

private:

	const double Lref_, Vref_, Tref_, qom_, Oref_;
	const IR3field *electric_field_, *magnetic_field_;
	const double iEfield_time_factor_, iBfield_time_factor_;
	const metric_connected *metric_;
	const morphism *my_morphism_;

}; // end class curvilinear_boris

} // end namespace gyronimo

#endif // GYRONIMO_CURVILINEAR_BORIS