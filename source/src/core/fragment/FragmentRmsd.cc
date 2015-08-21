// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/fragment/FragmentRmsd.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <core/fragment/FragmentRmsd.hh>

// Utility headers
#include <utility/exit.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

namespace core {
namespace fragment {

/// @details Auto-generated virtual destructor
FragmentRmsd::~FragmentRmsd() {}

FragmentRmsd::FragmentRmsd(FragSetCOP fragments) : fragments_(fragments) {
	for ( ConstFrameIterator i = fragments_->begin(); i != fragments_->end(); ++i ) {
		FrameCOP frame = *i;
		frames_[frame->start()] = frame;
	}
}

FrameCOP FragmentRmsd::frame(core::Size position) const {
	if ( frames_.find(position) == frames_.end() ) {
		utility_exit_with_message("Requested invalid position in FragmentRmsd::fragment");
	}

	return frames_[position];
}

FragDataCOP FragmentRmsd::fragment(core::Size position, core::Size k) const {
	FrameCOP f = frame(position);

	if ( k < 1 || k > f->nr_frags() ) {
		utility_exit_with_message("Requested invalid fragment number in FragmentRmsd::fragment");
	}

	return f->fragment_ptr(k);
}

core::Real FragmentRmsd::rmsd(core::Size position, core::Size k, const core::pose::Pose& reference) const {
	core::pose::Pose pose;
	core::pose::make_pose_from_sequence(pose, reference.sequence(), "centroid");

	FrameCOP f = frame(position);
	f->apply(k, pose);

	return core::scoring::CA_rmsd(pose, reference, f->start(), f->stop());
}

}  // namespace fragment
}  // namespace core
