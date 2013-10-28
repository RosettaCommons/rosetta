// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoopInserter.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_devel_loop_creation_LoopInserter_HH
#define INCLUDED_devel_loop_creation_LoopInserter_HH

//Unit
#include <devel/loop_creation/LoopInserter.fwd.hh>

//protocols
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>

//utility
#include <utility/tag/Tag.fwd.hh>

//C++
#include <set>

namespace devel {
namespace loop_creation {

class LoopInserter : public protocols::moves::Mover
{

public:
	
	protocols::loops::Loop
	get_created_loop() const;
	
	core::Size
	loop_anchor() const;
	
	void
	loop_anchor(
		core::Size loop_anchor
	);
	
	void
	modified_range(
		core::Size res_begin,
		core::Size res_end
	);
	
	std::pair<core::Size, core::Size>
	modified_range() const;
	
	void
	parse_loop_anchor(
		utility::tag::TagCOP const tag
	);
	
protected:

	//The inserted loop to be set by subclasses
	protocols::loops::Loop created_loop_;
	
	//Residues that have been modified by this loop inserter
	std::pair<core::Size, core::Size> modified_range_;
	
	//The residue to build the loop after (new loop will be between current loop_anchor_ and loop_anchor_+1)
	core::Size loop_anchor_;
	
	//Should this mover be allowed to make modifications outside of the loop region?
	bool prevent_nonloop_modifications_;
	
};

} //loop creation
} //devel

#endif
