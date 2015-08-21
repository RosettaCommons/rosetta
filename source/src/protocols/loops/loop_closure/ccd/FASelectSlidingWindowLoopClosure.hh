// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh
/// @brief header file for SlidingWindowLoopClosure protocol
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_FASelectSlidingWindowLoopClosure_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_FASelectSlidingWindowLoopClosure_hh

// Unit Headers
#include <protocols/loops/loop_closure/ccd/FASelectSlidingWindowLoopClosure.fwd.hh>
#include <protocols/loops/loop_closure/ccd/WidthFirstSlidingWindowLoopClosure.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/OrderedFragSet.hh>

//#include <protocols/simple_moves/FragmentMover.hh>
//#include <core/fragment/SecondaryStructure.hh>
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

/// Move these forward declarations to FASelectSlidingWindowLoopClosure.fwd.hh
class FASelectSlidingWindowLoopClosure;
typedef utility::pointer::shared_ptr< FASelectSlidingWindowLoopClosure > FASelectSlidingWindowLoopClosureOP;
typedef utility::pointer::shared_ptr< FASelectSlidingWindowLoopClosure const > FASelectSlidingWindowLoopClosureCOP;

class FASelectSlidingWindowLoopClosure : public WidthFirstSlidingWindowLoopClosure {
	typedef  WidthFirstSlidingWindowLoopClosure Parent;
public:
	/// @brief constructor: supply fragsets for fragment moves
	FASelectSlidingWindowLoopClosure(
		core::fragment::FragSetCOP fragset,
		core::scoring::ScoreFunctionOP scorefxn,
		core::kinematics::MoveMapCOP movemap
	);
	//@brief just set defaults -- expects fragset, scorefxn and movemap to be set later
	FASelectSlidingWindowLoopClosure();

	~FASelectSlidingWindowLoopClosure();
	virtual std::string get_name() const;

	static void register_options();
	void set_defaults();

	//@brief run find fragments that close loop  (if ideal loop closing: such that the less_cut pose is close RMSD <0.1 to pose more_cut)
	// returns less_cut and more_cut with best fragment already applied..
	virtual void select_final_loop( core::pose::Pose& more_cut, core::pose::Pose& less_cut );

	core::Real fascore( core::pose::Pose& fa_pose ) const;

	void set_fullatom_pose( core::pose::Pose& fa_pose );

private:
	core::pose::PoseOP fa_pose_;
};

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif  //INCLUDED_protocols_loops_loop_closure_ccd_FASelectSlidingWindowLoopClosure_hh
