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

// C/C++ headers
#include <iostream>

// Utility headers
#include <basic/Tracer.hh>

// Project headers
#include <core/types.hh>
//#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

namespace core {
namespace fragment {

static basic::Tracer TR("core.fragment.FragmentRmsd");

FragmentRmsd::FragmentRmsd(FragSetCOP fragments) : fragments_(fragments) {
  for (FrameIterator i = fragments_->begin(); i != fragments_->end(); ++i) {
    Frame const * frame = *i;
    frames_[frame->start()] = frame;
  }
}

core::Real FragmentRmsd::rmsd(core::Size position, core::Size k, const core::pose::Pose& reference) const {
  if (frames_.find(position) == frames_.end()) {  // invalid position
    TR.Warning << "Invalid position-- " << position << std::endl;
    return -1;
  }

  Frame const * frame = frames_[position];

  if (k < 1 || k > frame->nr_frags()) {  // invalid fragment
    TR.Warning << "Invalid fragment number-- " << k << std::endl;
    return -1;
  }

  core::pose::Pose pose;
  core::pose::make_pose_from_sequence(pose, reference.sequence(), "centroid");
  frame->apply(k, pose);

  //const FragData& fragment = frame->fragment(k);
  //fragment.apply(pose, frame->start(), frame->end());

  return core::scoring::CA_rmsd(pose, reference, frame->start(), frame->stop());
}

}  // namespace fragment
}  // namespace core
