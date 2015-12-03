// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/util.cc
/// @brief  Simple utilities for computing rotamer recovery between 2 poses
/// @author JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <core/pack/util.hh>

// Package Headers


// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <numeric/numeric.functions.hh>

// util
#include <core/pose/util.hh>
#include <utility/string_util.hh>
#include <core/types.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

namespace core {
namespace pack {

using namespace core;
using namespace core::pose;

static basic::Tracer TR( "core.pack.util" );

////////////////////////////////////////////////////////////////////////////////

/// @brief Percentage of residues which have same rotamers
core::Real residue_rotamer_recovery( Pose & pose, Pose & ref_pose, core::Real difference ) {

	using namespace numeric;

	// get angle differences
	utility::vector1< utility::vector1< Real > > rotamer_diffs = get_rotamer_angle_diffs( pose, ref_pose );

	// initialize variables
	core::Real recovered( 0.0 );
	core::Real nres( static_cast< Real >( nres_protein( pose ) ) );

	// go through residues
	for ( core::Size i = 1; i <= nres; ++i ) {

		core::Real identical( 0.0 );

		// go through rotamers
		for ( core::Size j = 1; j <= rotamer_diffs[ i ].size(); ++j ) {

			// if difference is < difference degrees: recovered
			if ( abs_difference( rotamer_diffs[ i ][ j ], 0.0 ) <= difference || abs_difference( rotamer_diffs[ i ][ j ], 360.0 ) <= difference ) {
				identical += 1.0;
			}
		}

		// if all rotamers for this residue are identical, then they are recovered
		if ( identical == rotamer_diffs[ i ].size() ) {
			recovered++;
		}
	}

	TR << "recovered: " << recovered << std::endl;
	TR << "nres: " << nres << std::endl;

	// compute percentage
	core::Real percent_recovered = recovered / nres;
	return percent_recovered;

} // rotamer recovery

////////////////////////////////////////////////////////////////////////////////

/// @brief Percentage of rotamers recovered
core::Real rotamer_recovery( Pose & pose, Pose & ref_pose, core::Real difference ) {

	using namespace numeric;

	// get angle differences
	utility::vector1< utility::vector1< Real > > rotamer_diffs = get_rotamer_angle_diffs( pose, ref_pose );

	// initialize variables
	core::Real nrotamers( 0.0 );
	core::Real recovered( 0.0 );
	core::Real nres( static_cast< core::Real > ( nres_protein( pose ) ) );

	// go through residues
	for ( core::Size i = 1; i <= nres; ++i ) {

		// go through rotamers
		for ( core::Size j = 1; j <= rotamer_diffs[ i ].size(); ++j ) {
			nrotamers += 1.0;

			// if difference is < 5 degrees: recovered
			if ( abs_difference( rotamer_diffs[ i ][ j ], 0.0 ) <= difference || abs_difference( rotamer_diffs[ i ][ j ], 360.0 ) <= difference ) {
				recovered += 1.0;
			}
		}
	}

	TR << "recovered: " << recovered << std::endl;
	TR << "nrotamers: " << nrotamers << std::endl;

	// compute percentage
	core::Real percent_recovered = recovered / nrotamers;
	return percent_recovered;

} // rotamer recovery

////////////////////////////////////////////////////////////////////////////////

/// @brief Get rotamer angle differences
/// @details Outer vector is pose length, inner vector is different chi's
utility::vector1< utility::vector1< Real > > get_rotamer_angle_diffs( Pose & pose, Pose & ref_pose ) {

	// pose length check, ignores virtuals (such as MEM)
	if ( nres_protein( pose ) != nres_protein( ref_pose ) ) {
		utility_exit_with_message( "Can't compare poses of different lengths. Quitting." );
	}

	// pose sequence check
	if ( pose.sequence() != ref_pose.sequence() ) {

		// get rid of the membrane residue in pose or native
		std::string trimmed_pose_seq( pose.sequence() );
		std::string trimmed_ref_seq( ref_pose.sequence() );
		if ( utility::endswith( pose.sequence(), "Z" ) ) {
			trimmed_pose_seq = utility::trim( pose.sequence(), "Z" );
		}
		if ( utility::endswith( ref_pose.sequence(), "Z" ) ) {
			trimmed_ref_seq = utility::trim( ref_pose.sequence(), "Z" );
		}

		// if sequences are still different, replace "Z" in the sequence
		if ( trimmed_pose_seq.size() != trimmed_ref_seq.size() ) {

			// replacing Z in the sequence
			TR << "WARNING: replacing Z in the sequence. Might be a membrane residue or an inserted or unrecognized residue!" << std::endl;

			std::string replaced_pose_seq = utility::replace_in( trimmed_pose_seq, "Z", "" );
			std::string replaced_ref_seq = utility::replace_in( trimmed_ref_seq, "Z", "" );

			// if sequences are still different, cry
			if ( replaced_pose_seq != replaced_ref_seq ) {
				TR << "pose.sequence " << replaced_pose_seq << std::endl;
				TR << "refpose.sequence " << replaced_ref_seq << std::endl;
				utility_exit_with_message( "Can't compare rotamers of different sequences. Quitting." );
			}
		}
	}

	// initialize empty vector
	utility::vector1< utility::vector1< Real > > rotamer_angle_diffs;

	// iterate through residues
	core::Size i = 0; // for pose
	core::Size j = 0; // for ref pose
	while ( i < pose.total_residue() && j < ref_pose.total_residue() ) {

		++i;;
		++j;

		// skip the membrane residue(s)
		if ( pose.residue( i ).name3() == "MEM" ) {
			++i;
		}
		if ( ref_pose.residue( j ).name3() == "MEM" ) {
			++j;
		}

		if ( pose.residue( i ).name3() != ref_pose.residue( j ).name3() ) {
			TR << "WARNING: COMPARING TWO DIFFERENT RESIDUES " << pose.residue( i ).name3() << i << " and " << ref_pose.residue( j ).name3() << j << std::endl;
		}

		// get rotamers
		utility::vector1< Real > pose_chi = pose.residue( i ).chi();
		utility::vector1< Real > ref_chi = ref_pose.residue( j ).chi();
		utility::vector1< Real > chi_diff;

		// iterate through rotamers
		for ( core::Size k = 1; k <= pose_chi.size(); ++k ) {
			chi_diff.push_back( pose_chi[ k ] - ref_chi[ k ] );
		}

		rotamer_angle_diffs.push_back( chi_diff );
	}

	return rotamer_angle_diffs;

} // get rotamer angle diffs


} // namespace pack
} // namespace core
