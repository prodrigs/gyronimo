
#ifndef GYRONIMO_MORPHISM_SPHERICAL
#define GYRONIMO_MORPHISM_SPHERICAL

#include <cmath>

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

class morphism_spherical : public morphism {

public:

	//! Class Constructor
	morphism_spherical() : morphism() {};

	//! Class Destructor
	~morphism_spherical() {};

	//! Maps the curvilinear coordinates `q` into cartesian coordinates `x`.
	IR3 operator()(const IR3 &q) const override;

	//! Inverse transform from cartesian coordinates `x` into curvilinear coordinates `q`.
	IR3 inverse(const IR3 &x) const override;

	//! Returns the morphism's derivatives, correspondent to the covariant basis vectors in point `q`.
	dIR3 del(const IR3 &q) const override;
	dSM3 g_del(const IR3 &q) const override;

	//! Returns the jacobian of the transformation in point `q`.
	double jacobian(const IR3 &q) const override;

	//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
	dIR3 del_inverse(const IR3 &q) const override;

private:

}; // end class morphism_spherical

} // end namespace gyronimo

#endif // GYRONIMO_MORPHISM_SPHERICAL