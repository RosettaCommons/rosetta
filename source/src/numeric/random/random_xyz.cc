// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/random/random_xyz.hh
/// @brief  Random vectors and stuff
/// @author Will Sheffler

#include <numeric/random/random_xyz.hh>

namespace numeric {
namespace random {

using numeric::Real;
using namespace numeric;
using namespace numeric::random;

xyzVector<Real> random_vector_spherical(){
	return xyzVector<Real>(gaussian(),gaussian(),gaussian());
}
xyzVector<Real> random_vector_unit_cube(){
	return xyzVector<Real>(uniform(),uniform(),uniform());

}
xyzVector<Real> random_vector(){
	return random_vector_spherical();
}
xyzVector<Real> random_normal(){
	return random_vector_spherical().normalized();
}

// cheaper to use gaussian, or do sines?
Quaternion<Real> random_unit_quaternion(){
	Real u1=uniform();
	Real u2=uniform() * constants::d::pi_2;
	Real u3=uniform() * constants::d::pi_2;
	Real sinu2 = sin(u2);
	Real sinu3 = sin(u3);
	Real cosu2 = sqrt(1.0-sinu2*sinu2);
	Real cosu3 = sqrt(1.0-sinu3*sinu3);
	return Quaternion<Real>( sqrt(1-u1)*sinu2,
		sqrt(1-u1)*cosu2,
		sqrt(  u1)*sinu3,
		sqrt(  u1)*cosu3 );
}
xyzMatrix<Real> random_rotation(){
	Quaternion<Real> q = random_unit_quaternion();
	return xyzMatrix<Real>::cols(
		1.0 - 2.0*q.y()*q.y() - 2.0*q.z()*q.z(),       2.0*q.x()*q.y() - 2.0*q.z()*q.w(),       2.0*q.x()*q.z() + 2.0*q.y()*q.w(),
		2.0*q.x()*q.y() + 2.0*q.z()*q.w(), 1.0 - 2.0*q.x()*q.x() - 2.0*q.z()*q.z(),       2.0*q.y()*q.z() - 2.0*q.x()*q.w(),
		2.0*q.x()*q.z() - 2.0*q.y()*q.w(),       2.0*q.y()*q.z() + 2.0*q.x()*q.w(), 1.0 - 2.0*q.x()*q.x() - 2.0*q.y()*q.y()
	);
}
xyzTransform<Real> random_xform(){
	return xyzTransform<Real>(random_rotation(),random_vector_spherical());
}

xyzTransform<Real> gaussian_random_xform(Real const & angsd, Real const & movsd){
	return xyzTransform<Real>(
		rotation_matrix_degrees( xyzVector<Real>(gaussian(),gaussian(),gaussian()), gaussian()*angsd ),
		movsd*xyzVector<Real>( gaussian(),gaussian(),gaussian()) );
}


} // namespace random
} // namespace numeric

