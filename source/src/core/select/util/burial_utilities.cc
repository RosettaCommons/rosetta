// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util/burial_utilities.cc
/// @brief Utilities for determining whether a point is within a cone.
/// @author Gabe Rocklin (sidechain neighbour selection)
/// @author Vikram K. Mulligan (vmullig@uw.edu -- moving this class to core and refactoring for noncanonicals; added determine_whether_point_is_buried() function.)

#include <core/select/util/burial_utilities.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>

namespace core {
namespace select {
namespace util {

#define NFOLD_REDUCTION 10 //A tenfold reduction in the distance function value is the point at which we no longer consider whether we might be in the distance cone.

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
) {
	core::Real accumulator(0.0);
	core::Real const dist_cutoff_sq( std::pow(std::log( NFOLD_REDUCTION - 1 ) + dist_exponent*dist_midpoint, 2) ); //Used to skip calculation for distant points.

	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		core::chemical::ResidueType const &restype( pose.residue_type(ir) );
		core::id::AtomID at2( restype.aa() == core::chemical::aa_gly ? restype.atom_index("2HA") : restype.first_sidechain_atom(), ir );
		core::id::AtomID at1( restype.icoor(at2.atomno()).stub_atom1().atomno(), ir );
		numeric::xyzVector< core::Real > const conevect_coordinate_2( pose.xyz(at2) );

		//Some crude distance checks to accelerate this:
		if ( point_coordinates.distance_squared( conevect_coordinate_2 ) > dist_cutoff_sq ) continue;

		//Do the actual calculation:
		numeric::xyzVector< core::Real > const conevect( conevect_coordinate_2 - pose.xyz(at1) );
		accumulator += calculate_point_in_cone( point_coordinates, conevect, conevect_coordinate_2, angle_exponent, angle_shift_factor, dist_exponent, dist_midpoint );
		if ( accumulator > burial_threshold ) return true;
	}
	return false;
}

} //namespace util
} //namespace select
} //namespace core
