// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util/burial_utilities.hh
/// @brief Utilities for determining whether a point is within a cone.
/// @author Gabe Rocklin (sidechain neighbour selection)
/// @author Vikram K. Mulligan (vmullig@uw.edu -- moving this class to core and refactoring for noncanonicals; added determine_whether_point_is_buried() function.)

#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <core/pose/Pose.fwd.hh>

#ifndef INCLUDED_core_select_util_burial_utilities_HH
#define INCLUDED_core_select_util_burial_utilities_HH

namespace core {
namespace select {
namespace util {

/// @brief Given a point in 3D space, and a vector and floats defining a cone, determine the extent to which the point
/// is in the cone.
/// @details The return value ranges from 0 (not in the cone) to 1 (fully in the cone).  The cone has fuzzy boundaries, so
/// non-integer return values are possible.
/// @param [in] point_coordinates The coordinates of a test point that may or may not be in the cone.  For the layer selector, this is the beta- or alpha-carbon
/// of another residue.
/// @param[in] conevect A vector defining the direction in which the test cone opens.  For the layer selector, this is the CA-CB vector of an amino acid.
/// @param[in] conevect_coordinate_2 The coordinate in space of the base of the cone.  For the layer selector, this is the CB atom.
/// @param[in] angle_exponent A value defining how rapidly the cone falls off in the angular direction.
/// @param[in] angle_shift_factor A value that shifts the angluar falloff.
/// @param[in] dist_exponent A value defining how rapidly the cone falls off with distance.
/// @param[in] dist_midpoint A value defining the midpoint of the distance falloff.
inline core::Real calculate_point_in_cone(
	numeric::xyzVector<core::Real> const &point_coordinates,
	numeric::xyzVector<core::Real> const &conevect,
	numeric::xyzVector<core::Real> const &conevect_coordinate_2,
	core::Real const angle_exponent,
	core::Real const angle_shift_factor,
	core::Real const dist_exponent,
	core::Real const dist_midpoint
) {
	numeric::xyzVector< Real > vect(point_coordinates - conevect_coordinate_2);
	core::Real const dist_term(1.0 / (1.0 + exp( dist_exponent*(vect.length() - dist_midpoint)  ))); // dist_term = 1/(1+exp( n*(d - m))); sigmoidal falloff with midpoint at m, with sharpness controlled by n
	core::Real angle_term( ( conevect.dot(vect.normalize()) + angle_shift_factor ) / (1 + angle_shift_factor ) );
	if ( angle_term < 0 ) {
		angle_term = 0.0;
	}
	return (dist_term * pow(angle_term, angle_exponent) );
}

/// @brief Given a point in 3D space, determine whether or not that point is buried by the method of sidechain neighbor cones.
/// @note A crude distance cutoff metric is also used to make this calculation more efficient.
/// @details Returns true for burial, false otherwise.
/// @param[in] point_coordinates The 3D coordinates of the point in space.
/// @param[in] pose The pose, for reference.
/// @param[in] angle_exponent A value defining how rapidly the cone falls off in the angular direction.
/// @param[in] angle_shift_factor A value that shifts the angluar falloff.
/// @param[in] dist_exponent A value defining how rapidly the cone falls off with distance.
/// @param[in] dist_midpoint A value defining the midpoint of the distance falloff.
/// @param[in] burial_threshold The cutoff for considering a point to be buried.  Roughly, this is the number of cones that the point must be inside
/// in order to consider it "buried".
bool determine_whether_point_is_buried (
	numeric::xyzVector<core::Real> const &point_coordinates,
	core::pose::Pose const &pose,
	core::Real const angle_exponent,
	core::Real const angle_shift_factor,
	core::Real const dist_exponent,
	core::Real const dist_midpoint,
	core::Real const burial_threshold
);

} //namespace util
} //namespace select
} //namespace core

#endif // INCLUDED_core_select_util_burial_utilities_HH
