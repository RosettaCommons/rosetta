// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/SlidingWindowLoopClosure.hh
/// @brief header file for SlidingWindowLoopClosure protocol
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_WidthFirstSlidingWindowLoopClosure_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_WidthFirstSlidingWindowLoopClosure_hh

// Unit Headers
#include <protocols/loops/loop_closure/ccd/WidthFirstSlidingWindowLoopClosure.fwd.hh>
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/types.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/OrderedFragSet.hh>

//#include <protocols/evaluation/PoseEvaluator.hh>
// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

/// @detail
// ths derivation from SlidingWindowLoopClosure changes the order of sampling... instead of all attempts with a certain window size first
/// it will switch between different windows. If nothing good is found after one sweep of window sizes we bail out... (assuming that it is fruitless)
/// if things look promising more sampling time is spend to find a good loop
class WidthFirstSlidingWindowLoopClosure : public SlidingWindowLoopClosure {
	typedef SlidingWindowLoopClosure Parent;
public:
	/// @brief constructor: supply fragsets for fragment moves
	WidthFirstSlidingWindowLoopClosure(
		core::fragment::FragSetCOP fragset,
		core::scoring::ScoreFunctionOP scorefxn,
		core::kinematics::MoveMapCOP movemap
	);

	//@brief just set defaults -- expects fragset, scorefxn and movemap to be set later
	WidthFirstSlidingWindowLoopClosure();


	//@brief run find fragments that close loop  (if ideal loop closing: such that the less_cut pose is close RMSD <0.1 to pose more_cut)
	// returns less_cut and more_cut with best fragment already applied..
	void sample_loops( core::pose::Pose& more_cut, core::pose::Pose& less_cut ) override;
	// virtual void select_final_loop( core::pose::Pose& more_cut, core::pose::Pose& less_cut );
	static void register_options();
	std::string get_name() const override;

protected:

	void set_defaults();


	core::Real window_acceptance_ratio_; //=0.1; fast_closure must have at least 0.1*nr_fragments succesfully closed loops
	core::Size nr_scored_sampling_passes_;
	//                                         to accept this window for further scored (slower) sampling
	core::Size min_fast_loops_;
	core::Size min_breakout_fast_loops_;

	core::Size nr_scored_fragments_;
	core::Size give_up_;
};

} // nameaspce ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_WidthFirstSlidingWindowLoopClosure_hh
