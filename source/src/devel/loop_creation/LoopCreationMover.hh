// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoopCreationMover.hh
///
/// @brief A Mover that uses loophash to find fragments that can
/// bridge a gap with minimal modifications to the original pose.
/// @author Tim Jacobs


#ifndef INCLUDED_devel_loop_creation_LoopCreationMover_HH
#define INCLUDED_devel_loop_creation_LoopCreationMover_HH

//Unit
#include <devel/loop_creation/LoopCreationMover.fwd.hh>

//Protocols
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/simple_moves/TaskAwareMinMover.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <devel/loop_creation/LoopCloser.fwd.hh>
#include <devel/loop_creation/LoopInserter.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//Core
#include <core/pose/Pose.fwd.hh>

//Utility
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace loop_creation {
	
class LoopCreationMover : public protocols::moves::Mover {

public:

	LoopCreationMover();
	
	LoopCreationMover(
		LoopInserterOP loop_inserter,
		LoopCloserOP loop_closer,
		core::Size attempts_per_anchor,
		bool refine,
		bool design_loops,
		bool include_neighbors,
		bool minimize_loops
	);

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	
	void
	init();
	
	std::string get_name() const;

	void
	apply(
		core::pose::Pose & pose
	);
	
	void
	copy_last_loop_to_new_anchor(
		core::pose::Pose & pose,
		core::Size new_anchor
	);
	
	/// @brief design and minimize the loop region
	core::Real
	refine_loop(
		core::pose::Pose & pose,
		protocols::loops::Loop loop
	);
	
	/// @brief return the most recently created loop
	protocols::loops::Loop
	get_last_created_loop() const;
	
	LoopInserterOP
	loop_inserter() const;
	
	LoopCloserOP
	loop_closer() const;

	void
	update_anchors(
		utility::vector1<core::Size> & loop_anchors,
		protocols::loops::Loop const & new_loop,
		core::Size index_of_new_loop
	);
	
//	LoopRefinerOP
//	loop_refiner() const;
	
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
	);
	
private:
	
	void
	detect_loop_anchors(
		core::pose::Pose const & pose
	);
	
	
private:
	//LoopInserter is a mover subclass responsible for inserting the actual loop residues
	LoopInserterOP loop_inserter_;
	
	//LoopCloser is a mover subclass responsible for closing the loop once loop residues have been inserted
	LoopCloserOP loop_closer_;
	
	//LoopCloser is a mover subclass responsible for refining the closed loop (design, minimization, etc)
	//	LoopRefinerOP loop_refiner_;

	//Loop filter
	protocols::filters::FilterOP loop_filter_;

	//If a particular insertion+closure step fails to produce a structure
	//that passes the loop filter, how many times should we retry?
	core::Size attempts_per_anchor_;
	
	//Should any refinement steps be done (master bool, without this, no packing/minimization at all)
	bool refine_;

	//Should the loops be designed before evaluation?
	bool design_loops_;
	
	//Should loop-neighbors be included in packing/design?
	bool include_neighbors_;

	//Should the loops be minimized before evaluation?
	bool minimize_loops_;

	//Filter using the loop analyzer mover
	bool filter_by_lam_;
	core::Real lam_score_cutoff_;
	
	//The most recently created loop
	protocols::loops::Loop last_created_loop_;
	
	//To close multiple loops. This overrides the loop anchor set by the loop inserter
	utility::vector1<core::Size> loop_anchors_;

	//A hacky way to specify that you want the loop to be copied to other anchors (for fake symmetry)
	core::Size asym_size_;

	//DEBUG
	bool dump_pdbs_;
	

	core::scoring::ScoreFunctionOP scorefxn_;
};

} //loop creation
} //devel

#endif
