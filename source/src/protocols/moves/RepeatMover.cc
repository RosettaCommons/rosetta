// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/RepeatMover.cc
/// @brief Method definitions for RepeatMover

// Unit Headers
#include <protocols/moves/RepeatMover.hh>

// Utility Headers
#include <utility/vector1.hh>

// Basic Headers
#include <basic/Tracer.hh>


namespace protocols {
namespace moves {

using namespace core;

// Empty constructor (nmoves=1)
RepeatMover::RepeatMover() : Mover(), nmoves_(1) {}

/// @details repeatmover = RepeatMover( mover_in , nmoves_in )
///
/// Mover    mover_in    /object defining what move to make
/// int      nmoves_in   /how many times to apply mover_in
RepeatMover::RepeatMover(
		MoverOP mover_in,
		int nmoves_in
) : Mover("RepeatMover"),
		mover_(mover_in),
		nmoves_(nmoves_in)
{}

// Copy constructor
RepeatMover::RepeatMover(RepeatMover const & object_to_copy) : Mover(object_to_copy),
		mover_(object_to_copy.mover_),
		nmoves_(object_to_copy.nmoves_)
{}

RepeatMover::~RepeatMover() {}

/// @details
/// Example(s):
///     repeatmover.apply(pose)
/// See Also:
///     MinMover
///     SequenceMover
///     ShearMover
///     SmallMover
///     TrialMover
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

protocols::moves::MoverOP
RepeatMover::clone() const
{
	return protocols::moves::MoverOP( new RepeatMover(*this) );
}

protocols::moves::MoverOP
RepeatMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RepeatMover() );
}

core::Size
RepeatMover::get_nmoves() const {
	if ( mover_ != 0 ) {
		return nmoves_;
	}
	else { return 0; }
}

std::string
RepeatMover::get_mover() const {
	if ( mover_ != 0 ) {
		return mover_->get_name();
	}
	else { return "none";}
}

std::ostream &operator<< (std::ostream &os, RepeatMover const &mover)
{

	os << "Mover name: " << mover.get_name() << ", Mover type: " << mover.get_type() << ", Mover current tag: " << mover.get_current_tag() << "\n" <<
			"Mover being repeated: " << mover.get_mover() << ", nmoves: " << mover.get_nmoves() << "\n";
	return os;
}

} // moves
} // protocols
