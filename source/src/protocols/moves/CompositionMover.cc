// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/CompositionMover.cc
/// @brief
/// @author

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/CompositionMover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

CompositionMover::CompositionMover()
: Mover("CompositionMover")
{}

void
CompositionMover::apply( core::pose::Pose & pose ) {
	for ( Size n = 1; n <= movers_.size(); n++ ) apply_mover_if_defined( pose, n );
}

std::string
CompositionMover::get_name() const {
	return "CompositionMover";
}

void CompositionMover::clear() {
	movers_.clear();
}

void CompositionMover::add_mover( MoverOP m ) {
	movers_.push_back( m );
}

utility::vector1< MoverOP > CompositionMover::get_movers() {
	return movers_;
}

void
CompositionMover::apply_mover_if_defined( core::pose::Pose & pose, Size const n ){
	if ( movers_[n] ) movers_[n]->apply( pose );
}

void
CompositionMover::apply( core::pose::Pose & pose, Size const i, Size const j ){
	runtime_assert( i > 0 );
	runtime_assert( j > 0 );
	runtime_assert( i <= movers_.size() );
	runtime_assert( j <= movers_.size() );
	if ( i <= j ) {
		for ( Size n = i; n <= j; n++ ) apply_mover_if_defined( pose, n );
	} else {
		for ( Size n = i; n >= j; n-- ) apply_mover_if_defined( pose, n );
	}
}

} // moves
} // protocols

