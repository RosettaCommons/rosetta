// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/CAlpha.hh
/// @brief  information about the c-alphas
/// @author Monica Berrondo


#ifndef INCLUDED_core_kinematics_CAlpha_hh
#define INCLUDED_core_kinematics_CAlpha_hh


// Package headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Rosetta Headers

// ObjexxFCL Headers

// // C++ Headers
#include <vector>


namespace core {
namespace kinematics {


/// CAlpha pair distance information
class CAlpha
{
public:
	// WARNING: this is really inefficient, fix later!!!!
	std::map < std::pair < Size, Size >, Real > squared_distance;

	/// caculate all C-alpha atom pairwise distance from a pose
	void
	update_distances_from_pose
	(
		pose::Pose const & pose
	)
	{
		Size const ca ( 2 );
		for ( Size i=1; i<=pose.total_residue(); ++i ) {
			PointPosition const & i_pos ( pose.residue( i ).atom( ca ).xyz() );
			for ( Size j=i+1; j<=pose.total_residue();++j ) {
				DistanceSquared const d_sq ( i_pos.distance_squared( pose.residue( j ).atom( ca ).xyz() ) );
				squared_distance[ std::make_pair(i,j) ] = d_sq;
			} ///< for loop over j
		} ///< for loop over i
	}

	/// retrieve distance for a pair of C-alpha atoms.
	Real
	get_squared_distance
	(
		Size const i,
		Size const j
	)
	{
		return squared_distance[ std::make_pair(i,j) ];
	}

};


} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_CAlpha_HH
