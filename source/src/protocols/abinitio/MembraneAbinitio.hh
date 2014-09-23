// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/MembraneAbinitio.hh
/// @brief header file for MembraneAbinitio protocol
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_MembraneAbinitio_hh
#define INCLUDED_protocols_abinitio_MembraneAbinitio_hh

// Unit Headers

// Package Headers
#include <protocols/simple_moves/FragmentMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/MembraneTopology.fwd.hh>
#include <core/scoring/MembranePotential.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.fwd.hh>

#include <protocols/abinitio/Protocol.hh>
#include <protocols/moves/TrialMover.hh>

#include <protocols/checkpoint/CheckPointer.hh>

//// C++ headers
#include <string>

#include <core/fragment/FragSet.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

/// Move these forward declarations to MembraneAbinitio.fwd.hh
class MembraneAbinitio;
typedef utility::pointer::shared_ptr< MembraneAbinitio > MembraneAbinitioOP;

//@ brief The Classic Abinitio protocol from rosetta++
/*!
@ detail
general usage:
MembraneAbinitio  abinitio;
abinitio.init( pose );
...
while(nstruct) {
	 abinitio.apply( pose );
}

call MembraneAbinitio::register_options() before core::init::init to add relevant options to the applications help

, with the following
stages, all of which uses a different ScoreFunction based on the cen_std.wts in minirosetta_database:

- Stage 1: large (usually 9mer) randomly selected fragment insertions, only VDW term turned on.
Uses score0.wts_patch and runs for either a maximum of 2000 cycles or until all moveable phi/psi values
have been changed.

- Stage 2: large randomly selected fragment insertions, more score terms turned on. Uses score1.wts_patch
and runs for 2000 cycles.

- Stage 3: uses large randomly selected fragment insertions, although the size of the fragment insertions
is tunable via the set_apply_large_frags( bool ) method. Alternates between score2.wts_patch and score5.wts_patch,
running tunable numbers of 2000-cycle iterations between the two scoring functions.

- Stage 4: uses small (usually 3mer) fragment insertions with the fragment selection based on the Gunn cost for
finding local fragment moves. Runs for 4000-cycles and uses score3.wts_patch.

The class implements the basic abinito approach as known from rosetta++. We tried to set this up, such that
behaviour of the protocol can be changed in many different ways ( see, e.g., FoldConstraints ). To be able to change the
behaviour of the protocol easily the class-apply function and methods called therein (e.g., prepare_XXX() / do_XXX_cycles() ) should
not directly change moves or trials. A reference to the currently used score-function should be obtained by
mc().score_function() ...

Behaviour can be changed in the following ways:

use non-classic FragmentMover  --> eg. not uniformly sampled fragments, but using some weighting
															 --> large and small moves doesn't have to be 3mers and 9mers... use other movers...
															 ---> or other fragets for the "convenience constructor"
use custom trial classes --> overload update_moves()

change sampling behaviour:
	 overload prepare_XXX() methods: these are called before the cycling for a certain stage begins
	 overload do_stageX_cycles() : the actual loops over trial-moves ...

change scoring functions:
	 overload set_default_scores()
	 weight-changes effective for all stages: set_score_weight()

*/


class MembraneAbinitio : public Protocol {
	typedef Protocol Parent;
public:
	enum StageID {
		ALL_STAGES = 0,
		STAGE_1,
		STAGE_2,
		STAGE_3a,
		STAGE_3b,
		STAGE_4
	};
	///@brief This constructor does not work -- Fix it before using it.
	// constructor: supply mover classes for Fragment Moves
	MembraneAbinitio(
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_small_top25,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small,
		int  /*dummy otherwise the two constructors are ambiguous */
	);

	///@brief constructor: supply fragsets for large and small fragment moves
	MembraneAbinitio(
					 core::fragment::FragSetCOP fragset_small,
					 core::fragment::FragSetCOP fragset_small_top25,
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap
	);

	//@brief setup moves, mc-object, scores
	//@details can't call this from constructor; virtual functions don't operate
	//until construction has completed.
	virtual
	void init( core::pose::Pose const& pose );

	//@brief MembraneAbinitio has virtual functions... use this to obtain a new instance
	virtual
	moves::MoverOP clone() const;
	//	MembraneAbinitioOP clone() const;

	//@brief run protocol on pose
	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	//@brief return FramgentMover for smooth_small fragment insertions (i.e., stage4 moves)
	simple_moves::FragmentMoverOP smooth_move_small();
	simple_moves::FragmentMoverOP brute_move_small_top25();

	//@brief return FragmentMover for small fragment insertions ( i.e., stage3/4 moves )
	simple_moves::FragmentMoverOP brute_move_small();

	//@brief return FragmentMover for large fragment insertions (i.e., stage1/2 moves )
	simple_moves::FragmentMoverOP brute_move_large();

	//@brief change the movemap ( is propagated to mover-objects )
	//@detail overload if your extension stores additional moves as member variables
	virtual void set_movemap ( core::kinematics::MoveMapCOP mm );

	//@brief returns current movemap
	core::kinematics::MoveMapCOP movemap();


	//@brief add a new transmembrane region to be inserted
	void add_spanning_region( core::pose::Pose & pose);
	void move_all_inserted( core::pose::Pose & pose);
	void print_debug( core::pose::Pose & pose);

	//@brief set new instances of FragmentMovers
	void set_moves (
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small
	);

	//@brief set new monte-carlo object
	void set_mc( moves::MonteCarloOP );

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarlo & mc() {
		return *mc_;
	}

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarlo const& mc() const {
		return *mc_;
	}

	//@brief override cycle setting ( init() -> sets it according to cmd-line options )
	virtual void set_cycles( core::Real increase_cycles = 1.0 );

	//@brief
	Size total_trials() const {
		return total_trials_;
	}

	//@brief override cycle setting for specific stage ( valid nr: 1,2,3,4)
	//	void set_cycles( Size stage, Size cycles );

	// eventually replace with a polymorphic checkpointer class
	// these functions are almost identical to the ones used for ClassicRelax
	virtual void clear_checkpoints(void);
	virtual bool recover_checkpoint( core::pose::Pose &pose, std::string const & id );
	virtual void checkpoint( core::pose::Pose &pose, std::string const & id );

	/* obsoleted by set_current_tag in Mover class
 /// @brief set tag for debugging output
		void set_output_tag( std::string tag ) {
		output_tag_ = tag;
	}

	std::string get_output_tag() {
		return output_tag_;
	}
	*/

	/// @brief for debugging, one wants to have access to the native pose.
	//	void set_native_pose( core::pose::Pose const & pose );

	//@brief set weight - effective for all scoring functions  stage == -1 --> set weight for all stages
	// mod -1 (ignored) ,  in stage3 mod = 1 --> 3a , mod = 2 --> 3b
	virtual void set_score_weight( core::scoring::ScoreType, core::Real setting, StageID stage = ALL_STAGES );

protected:
	//@brief called to notify about changes regarding movers... new movemap / new instances of FragmentMover
	virtual void update_moves();

	//@brief called by init() --- calls all set_default_XXX methods
	virtual void set_defaults( core::pose::Pose const& pose );

	//@brief read out cmd-line options
	virtual void set_default_options();

	//@brief register cmd-line options in option system ( call before core::init::init )
public:
	static void register_options();
protected:
	//@brief construct default monto-carlo object
	virtual void set_default_mc(
		core::pose::Pose const& pose,
		core::scoring::ScoreFunction const& scorefxn
	);

	//@brief assigns default score-functions to the 4 stages--> score0 / score1 / (score2/score5) / score3
	virtual void set_default_scores();

	//@brief currently used score function ( depends on stage )
	core::scoring::ScoreFunction const& current_scorefxn() const;

	//@brief set current scorefunction
	void current_scorefxn( core::scoring::ScoreFunction const& scorefxn );

	//@brief set individual weight of current scorefunction
	//@details NOTE: does not change the predefined weights, this information is lost each time the score is set
	//  weight will be overwritten with default ->at beginning of stages, each iteration of loop in stage3
	// we could change this behaviour by having a pointer to the current score and making sure that nobody can
	// change the score inside the mc-object (only const accessor )
	void set_current_weight( core::scoring::ScoreType type, core::Real setting );

	//@brief run cycles for different scoring_stages, return number of steps used
	virtual int do_stage1_cycles( core::pose::Pose &pose );
	virtual int do_stage2_cycles( core::pose::Pose &pose );
	virtual int do_stage3_cycles( core::pose::Pose &pose );
	virtual int do_stage3b_cycles( core::pose::Pose &pose );
	virtual int do_stage4_cycles( core::pose::Pose &pose );

	//@brief returns true if pose is < 3.0 A rms to last pose sent to this function
	//	bool convergence_check( core::pose::Pose const & pose );

	//@brief returns the Mover that is applied inside the stage3 double loop
	virtual moves::TrialMoverOP
	stage2_mover( core::pose::Pose &pose,  moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage3 double loop
	virtual moves::TrialMoverOP
	stage3_mover( core::pose::Pose & pose, int lct1, int lct2,  moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage4 loop
	virtual moves::TrialMoverOP
	stage4_mover( core::pose::Pose &pose, int kk, moves::TrialMoverOP trials_in );


	//@brief called by update_moves() creates the instances of TrialMover with the FragmentMoves
	virtual void set_trials();

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_large() {
		return ( apply_large_frags_ ? trial_large_ : trial_small_ );
	}

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_small() {
		return trial_small_;
	}

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_small_top25() {
		return trial_small_top25_;
	}

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_smooth() {
		return smooth_trial_small_;
	}

	// anything you want to have done before the stages ?
	//@brief prepare_stageX is called before do_stageX_cycles... overload to change status/scoring/conformation....
	virtual void prepare_stage1( core::pose::Pose &pose );
	virtual void prepare_stage2( core::pose::Pose &pose );
	virtual void prepare_stage3( core::pose::Pose &pose );
	virtual void prepare_stage4( core::pose::Pose &pose );

	//@brief called in each iteration of inner loop in stage3 before stage3_cycles_ of trials commence
	virtual void prepare_loop_in_stage3(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);

	//@brief called in each iteration of the loop in stage4 before the stage4_cycles_ of trials commence
	virtual void prepare_loop_in_stage4(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);

public:
	//@brief accessor for boolean flag: just_smooth_cycles
	inline
	bool just_smooth_cycles() const {
		return just_smooth_cycles_;
	}

	//@brief Accessor for number of stage1 cycles
	inline
	Size stage1_cycles() const {
		return stage1_cycles_;
	}

	//@brief Accessor for number of stage2
	inline
	Size stage2_cycles() const {
		return stage2_cycles_;
	}

	//@brief Accessor for number of stage3 cycles
	inline
	Size stage3_cycles() const {
		return stage3_cycles_;
	}

	//@brief Setter for number of stage4 cycles
	inline
	void set_stage4_cycles(Size stage4_cycles_new) {
		stage4_cycles_ = stage4_cycles_new;
	}


	//@brief Accessor for number of stage4 cycles
	inline
	Size stage4_cycles() const {
		return stage4_cycles_;
	}

	//@brief query this flag if you have time-intensive stuff and cut it short
	bool bQuickTest() const {
		return bQuickTest_;
	}

protected:
	void output_debug_structure( core::pose::Pose & pose, std::string prefix );

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarloOP mc_ptr() {
		return mc_;
	}


protected:
	//@brief  cycles_numbers  -- modified by set_cycles(), set_default_options()
	Size stage1_cycles_; //score0
	Size stage2_cycles_; //score1
	Size stage3_cycles_; //score25
	Size stage4_cycles_; //score3

public:

	void set_only_stage4( bool value ){ only_stage4_ = value; } //just_smooth_cycles_=value; }


private:
	Size get_tmh_inserted(core::pose::Pose const & pose) const;
	core::scoring::MembraneEmbed const & MembraneEmbed_from_pose( core::pose::Pose const & ) const;
	core::scoring::MembraneEmbed & nonconst_MembraneEmbed_from_pose( core::pose::Pose & ) const;
	core::scoring::MembraneTopology const & MembraneTopology_from_pose( core::pose::Pose const & ) const;
	core::scoring::MembraneTopology & nonconst_MembraneTopology_from_pose( core::pose::Pose & ) const;

	//@brief scoring functions -- modified by set_default_scores() / set_score_weight()
	core::scoring::ScoreFunctionOP score_stage1_;  //score0
	core::scoring::ScoreFunctionOP score_stage2_;  //score1
	core::scoring::ScoreFunctionOP score_stage3a_; //score2
	core::scoring::ScoreFunctionOP score_stage3b_; //score5
	core::scoring::ScoreFunctionOP score_stage4_;  //score3

	//@brief flags
	bool apply_large_frags_;   //above contig_cut2
	bool short_insert_region_; //below contig_cut3
	bool just_smooth_cycles_;
	bool only_stage4_;
	bool bQuickTest_;  // this flag might land in base-class


	//@brief a temperature
	core::Real temperature_;

	//@brief movemap --> which dofs can be moved during abinitio
	core::kinematics::MoveMapCOP movemap_;

	//@brief a MonteCarlo object -- set_default_mc() , access: mc()
	moves::MonteCarloOP mc_;

	// Large and small Fragments
	simple_moves::FragmentMoverOP brute_move_small_;
	simple_moves::FragmentMoverOP brute_move_small_top25_;
	simple_moves::FragmentMoverOP brute_move_large_;
	simple_moves::FragmentMoverOP smooth_move_small_;

	// TrialMovers
	moves::TrialMoverOP trial_large_;
	moves::TrialMoverOP trial_small_;
	moves::TrialMoverOP trial_small_top25_;
	moves::TrialMoverOP smooth_trial_small_;

	std::vector< std::string > checkpoint_ids_;

	Size total_trials_;

	checkpoint::CheckPointer checkpoint_;
};

} // abinitio
} // protocols

#endif
