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


//namespace protocols {
//namespace loophash {
namespace devel {
namespace loop_creation {
	
class LoophashLoopInserter : public devel::loop_creation::LoopInserter
{
public:

	typedef std::map< core::Size, std::vector<core::Size> > HashBuckets;

	LoophashLoopInserter();
	
	protocols::moves::MoverOP
	clone() const;
	
	protocols::moves::MoverOP
	fresh_instance() const;
		
	std::string
	get_name() const;
	
	virtual void
	apply(
		core::pose::Pose & pose
	);
	
	HashBuckets
	find_fragments(
		core::pose::Pose const & pose,
		core::Size lh_fragment_begin,
		core::Size lh_fragment_end
	);
	
	///@brief return all loophash fragments within the max_lh radius that also satisfy
	///torsion rms to the flanking regions
	HashBuckets
	find_fragments(
		core::pose::Pose const & pose,
		core::Size lh_fragment_begin,
		core::Size lh_fragment_end,
		core::Size min_fragment_size,
		core::Size max_fragment_size
	);

	///@brief get a random fragment length and fragment retrieval index from the given hash bucket.
	std::pair<core::Size,core::Size>
	get_random_fragment(
		HashBuckets hash_buckets
	);
	
	void
	init(
		core::pose::Pose & pose
	);
	
	///@brief build the specified loophash fragment into the pose. Return
	///the deviation from an ideal bond across the cut
	std::pair<core::Real, core::Real>
	build_loop(
		core::pose::Pose & pose,
		core::Size lh_fragment_begin,
		core::Size lh_fragment_end,
		core::Size lh_fragment_size,
		core::Size retrieve_index
	);
	
	virtual void
	parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
	);
	
protected:

	protocols::loophash::LoopHashLibraryOP lh_library_;
//	HashBuckets hash_buckets_;
	
	core::Real min_torsion_rms_;
	core::Real max_torsion_rms_;
	core::Real max_lh_radius_;
	
	core::Real max_closure_deviation_;
	
	utility::vector1<core::Size> loop_sizes_;
	
	//number of flanking residues for each segment that must match the
	//torsions of the input pose (within loophash min and max rms)
	core::Size num_flanking_residues_to_match_;
	
	bool modify_flanking_regions_;
	
	bool lh_initialized_;
};

} //loop creation
} //devel
//} //protocols
//} //loophash

#endif
