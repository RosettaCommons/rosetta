// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/NullDesignMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/movers/NullDesignMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

NullDesignMover::~NullDesignMover() {}

void
NullDesignMover::apply( core::pose::Pose & )
{
}

std::string
NullDesignMover::get_name() const {
	return "Null"; /// this actually conflicts with the other NullMover...
}

} //movers
} //protein_interface_design
} //protocols
