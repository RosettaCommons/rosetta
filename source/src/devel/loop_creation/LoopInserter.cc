// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoopInserter.cc
///
/// @brief
/// @author Tim Jacobs

//Unit
#include <devel/loop_creation/LoopInserter.hh>

//Utility
#include <utility/tag/Tag.hh>

//C++
#include <set>

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
LoopInserter::modified_range(
	core::Size res_begin,
	core::Size res_end
){
	modified_range_=std::make_pair(res_begin, res_end);
}

std::pair<core::Size, core::Size>
LoopInserter::modified_range() const{
	return modified_range_;
}

void
LoopInserter::loop_anchor(
	core::Size loop_anchor
){
	loop_anchor_=loop_anchor;
}

void
LoopInserter::parse_loop_anchor(
	utility::tag::TagCOP tag
){
	if ( tag->hasOption("loop_anchor") ) {
		loop_anchor_ =
			tag->getOption<core::Size>("loop_anchor");
	}
	// else{
	//  utility_exit_with_message("You must specify the loop anchor to use for loop insertion with the 'loop_anchor' tag");
	// }
}

} //loop creation
} //devel
