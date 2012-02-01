// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief
/// @author Javier Castellanos	(javiercv@uw.edu)

// Unit header
#include <devel/constrained_sequence_design/SequenceConstraint.hh>

// Package headers

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

typedef core::Real Real;
typedef core::Size Size;

namespace devel {
namespace constrained_sequence_design {

SequenceConstraint::SequenceConstraint():
	pose_(0)
{
}

SequenceConstraint::~SequenceConstraint()
{
}

SequenceConstraint::SequenceConstraint(PoseOP pose)
		: pose_( *pose ) {
}

Real
SequenceConstraint::raw_score() {
		Real total = 0.0;
		for(Size i = 1; i <= pose_->total_residue(); ++i) {
				total += this->apply(pose_->aa(i) , i);
		}
		return total;
} // raw score
	
Real
SequenceConstraint::score() {
		return this->raw_score() * weight_;
} // score

void
SequenceConstraint::set_pose( Pose& p ){
	pose_ = PoseAP(&p);
} // set_pose

void
SequenceConstraint:: parse_my_tag(
utility::tag::TagPtr const tag,
protocols::moves::DataMap & data,
Pose const & pose )
{
	//*pose_ = pose;
} // parse_my_tag

} // constrained_sequence_design
} // devel

