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

// @boris_push.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_BORIS_PUSH
#define GYRONIMO_BORIS_PUSH

#include <gyronimo/core/IR3algebra.hh>

namespace gyronimo {

//! Gyronimo implementation of the Cartesian Boris push.
/*!
	This class implements a single step in the Boris push [C. K. Birdsall and 
	A. B. Langdon, Plasma Physics via Computer Simulation, CRC Press, 1991], 
	defined only for cartesian coordinates.
	This method assumes the electric field is zero.
	The boris push refers to a particular way of discretizing the equation 
	of motion for the velocity, given by:
	@f$
		\frac{\tilde{\textbf{v}}(\tau+\Delta\tau/2)-
				\tilde{\textbf{v}}(\tau-\Delta\tau/2)}{\Delta\tau} 
			= \Omega_{ref} \frac{\tilde{\textbf{v}}(\tau+\Delta\tau/2)
				+\tilde{\textbf{v}}(\tau-\Delta\tau/2)}{2} 
				\times \tilde{\textbf{B}}(\tau)
    @f$
	The following variables are adimensional: the time is normalized to `Tref` 
	with @f$ t=T_{ref}\,\tau @f$, the velocity to `Vref`=`Lref`/`Tref` with 
	@f$ \textbf{v}=V_{ref}\,\tilde{\textbf{v}} @f$ and the magnetic field 
	to `Bref` with @f$ \tilde{\textbf{B}} = \textbf{B}/B_{ref} @f$.
	In this normalization, the reference frequency factor is given by
	@f$ \Omega_{ref} = \frac{q_s \, B_{ref}}{m_s} \, T_{ref} @f$.
*/
IR3 boris_push(const IR3& v_old, const double& Oref, 
		const double& Bmag, const IR3& Bversor, const double& dt);

//! Gyronimo implementation of the Cartesian Boris push.
/*!
	This class implements a single step in the Boris push [C. K. Birdsall and 
	A. B. Langdon, Plasma Physics via Computer Simulation, CRC Press, 1991], 
	defined only for cartesian coordinates.
	The boris push refers to a particular way of discretizing the equation 
	of motion for the velocity, given by:
	@f$
		\frac{\tilde{\textbf{v}}(\tau+\Delta\tau/2)-
				\tilde{\textbf{v}}(\tau-\Delta\tau/2)}{\Delta\tau} 
			= \Omega_{ref} \left( \tilde{\textbf{E}}(\tau) + 
				\frac{\tilde{\textbf{v}}(\tau+\Delta\tau/2)
				+\tilde{\textbf{v}}(\tau-\Delta\tau/2)}{2} 
				\times \tilde{\textbf{B}}(\tau) \right)
    @f$
	The following variables are adimensional: the time is normalized to `Tref` 
	with @f$ t=T_{ref}\,\tau @f$, the velocity to `Vref`=`Lref`/`Tref` with 
	@f$ \textbf{v}=V_{ref}\,\tilde{\textbf{v}} @f$, the electric field to 
	`Eref` with @f$ \tilde{\textbf{E}} = \textbf{E}/E_{ref} @f$ and the 
	magnetic field to `Bref` with @f$ \tilde{\textbf{B}} = \textbf{B}/B_{ref} @f$.
	In this normalization, the reference frequency factor is given by
	@f$ \Omega_{ref} = \frac{q_s \, B_{ref}}{m_s} \, T_{ref} @f$, and the 
	electric and magnetic fields are assumed to respect Faraday's law, with
	@f$ V_{ref} = E_{ref} / B_{ref} @f$.
*/
IR3 boris_push(const IR3& v_old, const double& Oref, const IR3& Efield, 
		const double& Bmag, const IR3& Bversor, const double& dt);

} // end namespace gyronimo

#endif // GYRONIMO_BORIS_PUSH