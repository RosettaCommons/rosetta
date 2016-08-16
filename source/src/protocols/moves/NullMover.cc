// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/NullMover.cc
/// @author Sarel Fleishman (sarelf@uw.edu)

// Unit Headers
#include <protocols/moves/NullMover.hh>

#include <utility/vector1.hh>


// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace moves {

NullMover::NullMover() :
	MoveMapMover( "NullMover" )
{
}

NullMover::~NullMover() {}

void
NullMover::apply( core::pose::Pose & )
{
}

std::string
NullMover::get_name() const {
	return "NullMover";
}

protocols::moves::MoverOP
NullMover::clone() const{
	return protocols::moves::MoverOP( new NullMover( *this ) );
}

protocols::moves::MoverOP
NullMover::fresh_instance() const{
	return protocols::moves::MoverOP( new NullMover );
}


} // moves
} // protocols

