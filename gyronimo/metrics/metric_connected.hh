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

// @metric_connected.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_CONNECTED
#define GYRONIMO_METRIC_CONNECTED

#include <gyronimo/core/error.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

//! Covariant metric that is connected to its respective `morphism`.
class metric_connected : public metric_covariant {

public:
	metric_connected(const morphism *morph);
	virtual ~metric_connected() override {};

	virtual SM3 operator()(const IR3& r) const override;
	virtual dSM3 del(const IR3& r) const override;
	virtual double jacobian(const IR3& r) const override;
	virtual IR3 del_jacobian(const IR3& r) const override;

	virtual ddIR3 christoffel_first_kind(const IR3& r) const override;
	virtual ddIR3 christoffel_second_kind(const IR3& r) const override;

	virtual const morphism* morph() const {return morph_;};

protected:
	const morphism *morph_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_METRIC_CONNECTED
