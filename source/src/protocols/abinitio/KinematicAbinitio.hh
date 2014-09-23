// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Abinitio - Type sampling controlled by KinematicControl object
/// @detailed
/// @author Oliver Lange
///


#ifndef INCLUDED_protocols_abinitio_KinematicAbinitio_hh
#define INCLUDED_protocols_abinitio_KinematicAbinitio_hh


// Unit Headers
#include <protocols/abinitio/KinematicAbinitio.fwd.hh>

#ifdef PYROSETTA
//#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#endif


// Package Headers
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/KinematicControl.hh>
#include <protocols/jumping/JumpSetup.hh>

#include <protocols/simple_moves/FragmentMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>



#include <core/fragment/FragSet.fwd.hh>

#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh>

#include <core/types.hh>

#include <core/scoring/constraints/ConstraintSet.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers



namespace protocols {
namespace abinitio {

class KinematicAbinitio : public FoldConstraints {
	typedef FoldConstraints Parent;
public:
	KinematicAbinitio(
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small,
		int dummy /* otherwise the two constructors are ambigous */
	);

	KinematicAbinitio(
		core::fragment::FragSetCOP fragset3mer,
		core::fragment::FragSetCOP fragset9mer,
		core::kinematics::MoveMapCOP movemap
	);

	~KinematicAbinitio();

	virtual moves::MoverOP clone() const {
		return moves::MoverOP( new KinematicAbinitio(*this) );
	}

	static void register_options();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual void set_max_seq_sep( core::pose::Pose& pose, Size setting );

	// @brief overload to do start extra-round of jump_cycles()
	virtual bool prepare_stage1( core::pose::Pose &pose );
	virtual bool prepare_stage2( core::pose::Pose &pose );
	virtual bool prepare_stage3( core::pose::Pose &pose );

	//@brief called in each iteration of inner loop in stage3 before stage3_cycles_ of trials commence
	virtual bool prepare_loop_in_stage3(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);

	//@brief called in each iteration of the loop in stage4 before the stage4_cycles_ of trials commence
	virtual bool prepare_loop_in_stage4(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);

	//@brief returns the Mover that is applied inside the stage3 double loop
	//overload from ClassicAbinitio to perform also some wobble-type random moves
	virtual moves::TrialMoverOP
	stage1_mover( core::pose::Pose &pose, moves::TrialMoverOP trials_in );


	//@brief returns the Mover that is applied inside the stage3 double loop
	//overload from ClassicAbinitio to perform also some wobble-type random moves
	virtual moves::TrialMoverOP
	stage2_mover( core::pose::Pose &pose, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage3 double loop
	//overload from ClassicAbinitio to perform also some wobble-type random moves
	virtual moves::TrialMoverOP
	stage3_mover( core::pose::Pose &pose, int lct1, int lct2, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage3 double loop
	//overload from ClassicAbinitio to perform also some wobble-type random moves
	virtual moves::TrialMoverOP
	stage4_mover( core::pose::Pose &pose, int kk, moves::TrialMoverOP trials_in );

	///@brief set the closure_protocol... if not set no closure...
	void
	closure_protocol( loops::loop_closure::ccd::SlidingWindowLoopClosureOP closure_protocol ) {
		closure_protocol_ = closure_protocol;
	}

protected:
	//@brief assigns default score-functions to the 4 stages--> score0 / score1 / (score2/score5) / score3
	virtual void set_default_scores();

	//@brief read out cmd-line options
	virtual void set_default_options();

	//@brief overload to set the kinematic control score-terms
	virtual void replace_scorefxn( core::pose::Pose& pose, StageID, core::Real intra_stage_progress );

private:
	//@brief helper method to create_bb_moves for the stage3_- and stage4_mover methods.
	moves::MoverOP
	create_bb_moves(
		core::pose::Pose &pose,
		moves::MoverOP std_moves,
		bool bLargeWobble,
		core::Real crank_up_angle
	);

	//@brief dump the "jump.log" files
	void dump_jump_log( core::pose::Pose& pose, std::string const& file_name );

	//@brief combines the jump_mover() with stdmove and returns the new Combimover
	moves::MoverOP
	create_jump_moves( moves::MoverOP stdmove );

	//@brief chainbreak weight is ramped up during stage3 and stage4
	bool bRampChainbreaks_;

	//@brief overlap chainbreak will be ramped in in stage4
	bool bOverlapChainbreaks_;

	//@brief if set we attempt loop-closing
	loops::loop_closure::ccd::SlidingWindowLoopClosureOP closure_protocol_;

	core::scoring::constraints::ConstraintSetOP full_constraint_set_;

};



class JumpingFoldConstraintsWrapper : public KinematicAbinitio {
	typedef KinematicAbinitio Parent;
public:
	JumpingFoldConstraintsWrapper(
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small,
		jumping::BaseJumpSetupOP jump_def,
		int dummy /* otherwise the two constructors are ambigous */
	);

	JumpingFoldConstraintsWrapper(
		core::fragment::FragSetCOP fragset3mer,
		core::fragment::FragSetCOP fragset9mer,
		core::kinematics::MoveMapCOP movemap,
		jumping::BaseJumpSetupOP jump_def
	);

	virtual moves::MoverOP clone() const {
		return moves::MoverOP( new JumpingFoldConstraintsWrapper(*this) );
	}

	//static void register_options();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

protected:

private:
	//@brief The generally possible jumps, --> used to generate JumpSamples for each run
	jumping::BaseJumpSetupOP jump_def_;
	//loops::SlidingWindowLoopClosureOP closure_protocol_;
};



} //abinitio
} // protocols

#endif

