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
/// @author Darwin Fu for uniform sphere sampling

#ifndef INCLUDED_numeric_sampling_random_xyz_hh
#define INCLUDED_numeric_sampling_random_xyz_hh


// Unit headers
#include <numeric/random/random.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/Quaternion.hh>
#include <numeric/trig.functions.hh>


namespace numeric {
namespace random {

/// @brief A random vector chosen uniformly from the ball (volume enclosed within a sphere) of the given radius around the origin.
numeric::xyzVector<numeric::Real> uniform_vector_sphere(numeric::Real radius = 1);

/// @brief A random vector chosen with spherical symmetry around the origin.
/// @details Actual distribution is a 3D gaussian with unit variance centered at the origin.
numeric::xyzVector<numeric::Real> random_vector();

/// @brief A random vector chosen with spherical symmetry around the origin.
/// @details Actual distribution is a 3D gaussian with unit variance centered at the origin.
numeric::xyzVector<numeric::Real> random_vector_spherical();

/// @brief A random vector chosen uniformly from within the volume of a unit cube with opposite verticies at (0,0,0) and (1,1,1)
numeric::xyzVector<numeric::Real> random_vector_unit_cube();

/// @brief A random vector chosens uniformly from the surface of a unit sphere centered on the origin.
numeric::xyzVector<numeric::Real> random_normal();

numeric::Quaternion<numeric::Real> random_unit_quaternion();
numeric::xyzMatrix<numeric::Real> random_rotation();
numeric::xyzTransform<numeric::Real> random_xform();
numeric::xyzTransform<numeric::Real> gaussian_random_xform(numeric::Real const & angsd, numeric::Real const & movsd);

} // namespace random
} // namespace numeric


#endif // INCLUDED_numeric_random_random_HH
