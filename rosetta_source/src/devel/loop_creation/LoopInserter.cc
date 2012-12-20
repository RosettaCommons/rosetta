// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoopInserter.cc
///
/// @brief
/// @author Tim Jacobs

//Unit
#include <devel/loop_creation/LoopInserter.hh>

namespace devel {
namespace loop_creation {

protocols::loops::Loop
LoopInserter::get_created_loop() const{
	return created_loop_;
}

core::Size
LoopInserter::loop_anchor() const{
	return loop_anchor_;
}

void
LoopInserter::loop_anchor(
	core::Size loop_anchor
){
	loop_anchor_=loop_anchor;
}
	
} //loop creation
} //devel
