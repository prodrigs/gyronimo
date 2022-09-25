// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

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

// @metric_nexus.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_NEXUS
#define GYRONIMO_METRIC_NEXUS

#include <gyronimo/core/error.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/metrics/metric_covariant.hh>

namespace gyronimo {

// Naming options: 
// metric_nexus 	// because of the connection between metric and morphism
// metric_generator // because it is generated from morphism
// metric_creator	// because it is created from morphism
// metric_founder	// because it is founded from morphism
// metric_calculator	// because it is calculated from morphism
// metric_accessory	// because it is an accessory of morphism

//! Covariant metric that is connected to its respective `morphism`.
class metric_nexus : public metric_covariant {

public:
	metric_nexus(const morphism *morph);
	virtual ~metric_nexus() override {};

	virtual SM3 operator()(const IR3& r) const override;
	virtual dSM3 del(const IR3& r) const override;
	virtual double jacobian(const IR3& r) const override;
	virtual IR3 del_jacobian(const IR3& r) const override;

	virtual ddIR3 christoffel_first_kind(const IR3& r) const override;
	virtual ddIR3 christoffel_second_kind(const IR3& r) const override;

	virtual const morphism* morph() const {return morph_;};

private:
	const morphism *morph_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_METRIC_NEXUS
