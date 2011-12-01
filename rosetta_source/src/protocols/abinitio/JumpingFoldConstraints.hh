// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FoldConstraints.hh
/// @brief Abinitio-Folding under (distance-)constraints
/// @detailed
/// extension of classic Foldconstraints Protocol to enable jumping
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_JumpingFoldConstraints_hh
#define INCLUDED_protocols_abinitio_JumpingFoldConstraints_hh


// Unit Headers
// #include <core/fragment/JumpingFoldConstraints.fwd.hh>

// Package Headers
#include <protocols/abinitio/FoldConstraints.hh>
// AUTO-REMOVED #include <protocols/abinitio/KinematicAbinitio.hh>
#include <protocols/jumping/JumpSetup.hh>

#include <protocols/basic_moves/FragmentMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>

#include <core/fragment/FragSet.fwd.hh>

#include <core/types.hh>

//Auto Headers
#include <core/fragment/FragSet.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers



namespace protocols {
namespace abinitio {


class JumpingFoldConstraints : public FoldConstraints {
  typedef FoldConstraints Parent;
public:
  JumpingFoldConstraints(
			 basic_moves::FragmentMoverOP brute_move_small,
			 basic_moves::FragmentMoverOP brute_move_large,
			 basic_moves::FragmentMoverOP smooth_move_small,
			 jumping::BaseJumpSetupOP jump_def,
			 int dummy /* otherwise the two constructors are ambigous */
  );

  JumpingFoldConstraints(
			 core::fragment::FragSetCOP fragset3mer,
			 core::fragment::FragSetCOP fragset9mer,
			 core::kinematics::MoveMapCOP movemap,
			 jumping::BaseJumpSetupOP jump_def
  );
  ~JumpingFoldConstraints() {};
  virtual moves::MoverOP clone() const {
    return new JumpingFoldConstraints(*this);
  }

	static void register_options();

  virtual void apply( core::pose::Pose & pose );
  virtual void set_max_seq_sep( core::pose::Pose& pose, Size setting );
  void set_native_pose( core::pose::Pose const& native_pose); //keeps a copy
  void set_defeat_purpose( bool setting ) {
    bChainbreaksAlwaysActive_ = setting;
  }

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


	//@brief cycles-trial moves of jump_frags_
	void jump_cycles( core::pose::Pose& pose, Size cycles );

	//@brief returns the Mover that is applied inside the stage3 double loop
	//overload from ClassicAbinitio to perform also some wobble-type random moves
	moves::TrialMoverOP
	stage2_mover( core::pose::Pose &pose, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage3 double loop
	//overload from ClassicAbinitio to perform also some wobble-type random moves
	moves::TrialMoverOP
	stage3_mover( core::pose::Pose &pose, int lct1, int lct2, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage3 double loop
	//overload from ClassicAbinitio to perform also some wobble-type random moves
	moves::TrialMoverOP
	stage4_mover( core::pose::Pose &pose, int kk, moves::TrialMoverOP trials_in );

protected:
	virtual void setup_default_min_move();

	//@brief assigns default score-functions to the 4 stages--> score0 / score1 / (score2/score5) / score3
	virtual void set_default_scores();

	//@brief read out cmd-line options
	virtual void set_default_options();

private:
	//@brief helper method to create_bb_moves for the stage3_- and stage4_mover methods.
	moves::MoverOP
	create_bb_moves(
  	core::pose::Pose &pose,
		moves::MoverOP std_moves,
		bool bLargeWobble,
		core::Real crank_up_angle
	);

	//@brief helper method to create_bb_moves for the stage3_- and stage4_mover methods.
	moves::MoverOP
	create_jump_moves( moves::MoverOP stdmove );

  void setup_foldtree( core::pose::Pose& );

  //@brief The generally possible jumps, --> used to generate JumpSamples for each run
	jumping::BaseJumpSetupOP jump_def_;

  //@brief if true chainbreaks are active from step 1
  bool bChainbreaksAlwaysActive_;

	//@brief chainbreak weight is ramped up during stage3 and stage4
	bool bRampChainbreaks_;

	//@brief overlap chainbreak will be ramped in in stage4
	bool bOverlapChainbreaks_;

	//@brief use jumps from ss-pair library
	bool bSampleJumps_;

	//@brief steal jump-geometry from native structure
	//( i.e., add them temporarily to ss-pair library )
	bool bStealJumps_;

	//@brief sample jumps from ss-library also during fragment insertion
	bool bJumpFragMoves_;

  //@brief the current set of jumps and chainbreaks
	jumping::JumpSample current_jumps_;

  core::pose::PoseOP native_pose_; // this is to be able to copy jumps from native structure
  core::pose::Pose starting_pose_; // this is to be able to rewind folding after jumps have been sampled

	core::fragment::FragSetOP pure_large_frags_; //keep these so that we can create
	core::fragment::FragSetOP jump_frags_; //keep these so that we can create

	Size nr_jump_cycles_; // controls number of pure-jump moves at beginning of each stage

};



} //abinitio
} // protocols

#endif


