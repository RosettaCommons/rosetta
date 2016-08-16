// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/import_pose/pose_stream/MetaPoseInputStream.cc
/// @brief
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {
namespace pose_stream {

void
MetaPoseInputStream::add_pose_input_stream( PoseInputStreamOP input ) {
	input_streams_.push_back( input );
}

utility::vector1< PoseInputStreamOP >
MetaPoseInputStream::get_input_streams() {
	return input_streams_;
}

bool MetaPoseInputStream::has_another_pose() {

	if ( current_index_ > input_streams_.size() ) return false;
	if ( input_streams_[ current_index_ ]->has_another_pose() ) return true;

	// move to next input_stream that has an available pose if necessary
	while ( current_index_ <= input_streams_.size() &&
			!input_streams_[current_index_]->has_another_pose()
			)
			++current_index_;

	if ( current_index_ > input_streams_.size() ) return false;

	// if we've gotten here, current_input_stream_ has another pose!
	return true;
}

void MetaPoseInputStream::reset(){
	for ( Size n = 1;  n <= input_streams_.size(); n++ ) input_streams_[n]->reset();
	current_index_ = 1;
}


void MetaPoseInputStream::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) {
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message( "MetaPoseInputStream: called fill_pose, but I have no more Poses!" );
	}

	// (*current_input_stream_)->fill_pose( pose, residue_set );
	input_streams_[ current_index_ ]->fill_pose( pose, residue_set );
}

void MetaPoseInputStream::fill_pose(
	core::pose::Pose & pose
) {
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message( "MetaPoseInputStream: called fill_pose, but I have no more Poses!" );
	}

	// (*current_input_stream_)->fill_pose( pose, residue_set );
	input_streams_[ current_index_ ]->fill_pose( pose );
}

} // pose_stream
} // import_pose
} // core
