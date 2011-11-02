// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/RepeatMover.cc
/// @brief
/// @author

// Unit Headers
#include <protocols/moves/RepeatMover.hh>

#include <utility/vector1.hh>


// Package headers

// Project headers

// tracer
// AUTO-REMOVED #include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace moves {

/// RepeatMover
RepeatMover::RepeatMover() : Mover(), nmoves_(1) {}

RepeatMover::RepeatMover(
	MoverOP mover_in,
	int nmoves_in
) : Mover("RepeatMover"),
		mover_(mover_in),
		nmoves_(nmoves_in)
{}

RepeatMover::~RepeatMover() {}

void
RepeatMover::apply( core::pose::Pose & pose ) {
	for ( int i=1; i<=nmoves_; ++i ) {
//		T("protocols.moves.RepeatMover") << "Move: " << i << "/" << nmoves_ << std::endl;
		mover_->apply( pose );
	}
}

std::string
RepeatMover::get_name() const {
	return "RepeatMover";
}

} // moves
} // protocols

