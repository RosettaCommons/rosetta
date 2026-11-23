// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/AdditionalPosesMover.hh
/// @brief This mover is an abstract base class for movers which
///   accept multiple poses from a previous mover.
///   (This is for integration with other processing machinery.)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/moves/AdditionalPosesMover.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace moves {

static basic::Tracer TR("protocols.moves.AdditionalPoseMover");

void
AdditionalPosesMover::set_previous_mover( protocols::moves::MoverOP const m ) {
	previous_mover_ = m;
}

bool
AdditionalPosesMover::previous_mover_set() const {
	return previous_mover_ != nullptr;
}

std::string
AdditionalPosesMover::previous_mover_name() const {
	if ( previous_mover_ ) {
		return previous_mover_->get_name();
	} else {
		return "NO PREVIOUS MOVER"; // As we're probably experiencing a bug.
	}
}

core::pose::PoseOP
AdditionalPosesMover::get_next_pose() const {
	if ( previous_mover_ ) {
		return previous_mover_->get_additional_output();
	} else {
		TR.Warning << "Attempting to get the next pose from a previous mover which doesn't exist." << std::endl;
		return nullptr;
	}
}

utility::vector1< core::pose::PoseOP >
AdditionalPosesMover::get_next_pose_batch(core::Size N) const {
	utility::vector1< core::pose::PoseOP > pose_batch;
	if ( previous_mover_ ) {
		while ( N == 0 || pose_batch.size() < N ) {
			core::pose::PoseOP next_pose = previous_mover_->get_additional_output();
			if ( next_pose == nullptr ) { break; }
			pose_batch.push_back(next_pose);
		}
	} else {
		TR.Warning << "Attempting to get the next pose batch from a previous mover which doesn't exist." << std::endl;
	}
	return pose_batch;
}


} //moves
} //protocols

