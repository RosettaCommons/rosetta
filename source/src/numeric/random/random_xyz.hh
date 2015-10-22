// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random/random_xyz.hh
/// @brief  Random vectors and stuff
/// @author Will Sheffler

#ifndef INCLUDED_numeric_sampling_random_xyz_hh
#define INCLUDED_numeric_sampling_random_xyz_hh


// Unit headers
#include <numeric/random/random.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/Quaternion.hh>


namespace numeric {
namespace random {

numeric::xyzVector<numeric::Real> random_vector();
numeric::xyzVector<numeric::Real> random_vector_spherical();
numeric::xyzVector<numeric::Real> random_vector_unit_cube();
numeric::xyzVector<numeric::Real> random_normal();
numeric::Quaternion<numeric::Real> random_unit_quaternion();
numeric::xyzMatrix<numeric::Real> random_rotation();
numeric::xyzTransform<numeric::Real> random_xform();
numeric::xyzTransform<numeric::Real> gaussian_random_xform(numeric::Real const & angsd, numeric::Real const & movsd);

} // namespace random
} // namespace numeric


#endif // INCLUDED_numeric_random_random_HH
