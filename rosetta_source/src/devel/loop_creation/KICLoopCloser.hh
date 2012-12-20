// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file KICLoopCloser.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_devel_loop_creation_KICLoopCloser_HH
#define INCLUDED_devel_loop_creation_KICLoopCloser_HH

//Unit
#include <devel/loop_creation/LoopCloser.hh>

//protocols
#include <protocols/loops/Loop.hh>

namespace devel {
namespace loop_creation {

class KICLoopCloser : public LoopCloser
{
public:

	KICLoopCloser();

	bool
	close_loop(
		core::pose::Pose & pose,
		protocols::loops::Loop loop
	);
	
	///@brief parse tag for use in RosettaScripts
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
	);
	
private:
	bool prevent_nonloop_changes_;
};

} //loop_creation
} //devel

#endif
