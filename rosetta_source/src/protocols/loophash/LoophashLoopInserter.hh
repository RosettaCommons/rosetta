// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoophashLoophashLoopInserter.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_loophash_LoophashLoopInserter_HH
#define INCLUDED_protocols_loophash_LoophashLoopInserter_HH

//Unit
#include <devel/loop_creation/LoopInserter.hh>


//Protocols
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loops/Loop.fwd.hh>

namespace protocols {
namespace loophash {

class LoophashLoopInserter : public protocols::loops::loop_creation::LoopInserter
{
public:

	LoophashLoopInserter();
	
	protocols::loops::Loop
	insert_loop(
		core::pose::Pose & pose,
		core::Size loop_anchor
	);
	
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);
	
private:

	LoopHashLibraryOP lh_library_;
	
	core::Real min_torsion_rms_;
	core::Real max_torsion_rms_;
	core::Real max_lh_radius_;
	
	utility::vector1<core::Size> loop_sizes_;
	
	//number of flanking residues for each segment that must match the
	//torsions of the input pose (within loophash min and max rms)
	core::Size num_flanking_residues_to_match_;
	
	bool modify_flanking_regions_;
};

} //protocols
} //loophash

#endif
