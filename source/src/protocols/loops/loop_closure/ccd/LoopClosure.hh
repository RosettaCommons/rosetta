// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/LoopClosure.hh
/// @brief header file for LoopClosure protocol
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_LoopClosure_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_LoopClosure_hh

// Unit Headers
#include <protocols/loops/loop_closure/ccd/LoopClosure.fwd.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/types.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.hh>

#include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers

#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

class LoopClosure : public utility::pointer::ReferenceCount {
public:
	/// @brief constructor: supply fragsets for fragment moves
	LoopClosure(
		core::fragment::FragSetCOP fragset,
		core::scoring::ScoreFunctionOP scorefxn,
		Loop loop_def,
		core::kinematics::MoveMapCOP movemap
	);

	//destructor
	virtual ~LoopClosure();

	//@brief run protocol on pose
	virtual bool apply( core::pose::Pose const& pose );

	//@brief return the list of collected fragments
	// fo the basic LoopClosure class this will contain only 1 Frame. could have returned the frame
	// but maybe it is worth to keep the interface more general ?
	core::fragment::FrameOP
	closure_fragments() const {
		return closure_frame_;
	}

	//@brief returns current movemap
	core::kinematics::MoveMapCOP movemap() const;

	//  //@brief set new monte-carlo object
	//   void set_mc( moves::MonteCarloOP mc ) {
	//     mc_ = mc;
	//   }

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarlo & mc() {
		return *mc_;
	}

	//   //@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	//   moves::MonteCarlo const& mc() const {
	//     return *mc_;
	//   }

	core::scoring::ScoreFunction const& scorefxn() {
		return *scorefxn_;
	}


	//@brief override cycle setting, sets nr_fragments to 100*ratio
	//  and trials to 20*loopsize*ratio
	void set_cycles( core::Real cycle_ratio = 1.0 );

	void set_nr_fragments( core::Size nr_fragments = 100 );

	core::Size nr_fragments() const {
		return nr_fragments_;
	}

	void ramp_chainbreak( bool setting = true ) {
		bRampChainbreak_ = setting;
	}

protected:
	//protected c'stor: set up relevant stuff
	// and call set_defaults from derived-class c'stor
	LoopClosure();

	void init();

	//   //@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	//   moves::MonteCarloOP mc_ptr() {
	//     return mc_;
	//   }

	//@brief inner-loop of fragment and ccd-moves
	virtual void do_frag_cycles( core::pose::Pose &pose ) const;

	/// @brief save the loop-fragment in closure_frames_
	virtual void catch_fragment( core::pose::Pose const& pose );

	//   core::fragment::FrameOP&
	//   closure_fragments() {
	//     return closure_frame_;
	//   }

	/// @brief replace scorefxn
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn );

	void set_loop( Loop const& loop_in ) {
		loop_ = loop_in;
	}

	Loop const& loop() const {
		return loop_;
	}

	void init_mc();

	void set_movemap( core::kinematics::MoveMapCOP mm );

	void set_fragset( core::fragment::FragSetCOP frags );

	//@brief override temperature setting
	void set_temperature( core::Real temperature ) {
		temperature_ = temperature;
	}

	virtual void ramp_chainbreak( core::Size iter, core::Size total ) const;

protected:
	Loop loop_;

	core::scoring::ScoreFunctionOP scorefxn_;  //score3

	//@brief a temperature
	core::Real temperature_;

	//@brief movemap --> which dofs can be moved during loops
	core::kinematics::MoveMapCOP movemap_;

	//@brief a MonteCarlo object -- set_default_mc() , access: mc()
	moves::MonteCarloOP mc_;

	core::fragment::FrameOP closure_frame_;

	core::Size nr_fragments_; //outer_cycles;
	core::Size cycles_; // trials per fragments

	simple_moves::FragmentMoverOP frag_mover_;
	moves::MoverOP ccd_mover_;

	core::fragment::FragSetCOP fragset_;

	bool bEnableCcdMoves_;
	bool bRampChainbreak_;
	core::Real final_weight_linear_chainbreak_;
	core::Real final_weight_overlap_chainbreak_;
};

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_LoopClosure_hh
