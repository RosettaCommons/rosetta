// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mover.cc
/// @brief Method code and full headers for Mover--
/// keeps heavily-included Mover.hh small and concise to maximize compiling
/// efficiency and to make the class definitions easier to read.
/// @author

// Unit Headers
#include <protocols/moves/WhileMover.hh>

// Package headers

// Project headers

// tracer
#include <basic/Tracer.hh>

#include <utility>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace moves {

/// WhileMover
PoseCondition::~PoseCondition() = default;

WhileMover::WhileMover() : Mover(), nmoves_(1) {}

WhileMover::WhileMover(
	MoverOP mover_in,
	core::Size nmoves_in,
	PoseConditionOP condition
) :
	Mover("ConditionalRepeatMover"),
	mover_(std::move(mover_in)),
	nmoves_(nmoves_in),
	p_cond_(std::move( condition ))
{}

WhileMover::~WhileMover() = default;

void
WhileMover::apply( core::pose::Pose & pose ) {
	PoseCondition &cond_( *p_cond_ );
	for ( Size i=0; ( i<nmoves_ ) && cond_( pose ); ++i ) {
		mover_->apply( pose );
	}
}

std::string
WhileMover::get_name() const {
	return "WhileMover";
}


} // moves
} // protocols

