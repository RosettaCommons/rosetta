// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CCDLoopCloser.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_devel_loop_creation_CCDLoopCloser_HH
#define INCLUDED_devel_loop_creation_CCDLoopCloser_HH

//Unit
#include <devel/loop_creation/LoopCloser.hh>

//core
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>

//protocols
#include <protocols/loops/Loop.hh>

namespace devel {
namespace loop_creation {

class CCDLoopCloser : public LoopCloser
{
public:

	///@brief default constructor
	CCDLoopCloser();
	
	///@brief explicit constructor
	CCDLoopCloser(
		core::Size max_closure_attempts,
		bool prevent_nonloop_modifications,
		core::Size max_ccd_moves_per_closure_attempt,
		core::Real tolerance,
		core::Real max_rama_score_increase,
		core::Real max_total_delta_helix,
		core::Real max_total_delta_strand,
		core::Real max_total_delta_loop,
		core::Real early_exit_cutoff
	);
	
	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;
		
	std::string
	get_name() const;
	
	void
	init();
	
	void
	apply ( core::pose::Pose & );
	
	///@brief parse tag for use in RosettaScripts
	void
	parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
	);
	
	void
	prepare_fold_tree(
		core::pose::Pose & pose
	);
	
private:

	//Max number of CCD calls to make before giving up
	core::Size max_closure_attempts_;
	
	//Should we set up a fold tree that prevents modification outside the loop region?
	bool prevent_nonloop_modifications_;

	//Number of CCD moves per call to CCD closure
	core::Size max_ccd_moves_per_closure_attempt_;

	//The threshold for whether or not we consider this loop closed
	core::Real tolerance_;

	//Inputs to ccd
	core::Real max_rama_score_increase_;
	core::Real max_total_delta_helix_;
	core::Real max_total_delta_strand_;
	core::Real max_total_delta_loop_;
	core::Real early_exit_cutoff_;
};

} //loop creation
} //devel

#endif
