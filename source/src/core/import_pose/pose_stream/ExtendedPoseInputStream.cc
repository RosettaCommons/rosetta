// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.hh>

#include <core/import_pose/pose_stream/ExtendedPoseInputStream.hh>

// C++ headers

#include <utility/exit.hh>

#include <core/pose/annotated_sequence.hh>


namespace core {
namespace import_pose {
namespace pose_stream {

bool ExtendedPoseInputStream::has_another_pose() {
	return ( current_n_ <= ntimes_ );
}

void ExtendedPoseInputStream::reset() {
	current_n_ = 1;
}

void ExtendedPoseInputStream::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set,
	bool const metapatches /*= true*/
) {
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message(
			"ExtendedPoseInputStream: called fill_pose, but I have no more Poses!"
		);
	}

	core::pose::make_pose_from_sequence(
		pose,
		seq_,
		residue_set,
		metapatches
	);

	for ( Size pos = 1; pos <= pose.size(); pos++ ) {
		pose.set_phi  ( pos, -150 );
		pose.set_psi  ( pos,  150 );
		pose.set_omega( pos,  180 );
	}

	++current_n_;
} // fill_pose

void ExtendedPoseInputStream::fill_pose(
	core::pose::Pose &,
	bool const //metapatches /*= true*/
) {
	utility_exit_with_message(
		"ExtendedPoseInputStream: called fill_pose, but without ResidueType Set"
	);
}

/// @brief Get a string describing the last pose and where it came from.
/// @details Input tag + filename.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
ExtendedPoseInputStream::get_last_pose_descriptor_string() const {
	runtime_assert_string_msg( current_n_ > 1, "Error in ExtendedPoseInputStream::get_last_pose_descriptor_string(): The fill_pose() function must be called at least once before calling this function." );
	char curindex_str[64];
	sprintf( curindex_str, "%04lu", current_n_ - 1 );
	std::string const outstr( "extended_pose_" + std::string(curindex_str) );
	return outstr;
}

} // pose_stream
} // import_pose
} // core
