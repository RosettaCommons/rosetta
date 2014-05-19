// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/FragmentSampler.hh
/// @brief header file for FragmentSampler protocol
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka


#ifndef INCLUDED_protocols_abinitio_FragmentSampler_hh
#define INCLUDED_protocols_abinitio_FragmentSampler_hh

#include <protocols/abinitio/FragmentSampler.fwd.hh>


// Unit Headers

//#include <core/io/silent/SilentFileData.hh>
//#include <core/io/silent/SilentFileData.fwd.hh>
//#include <core/io/silent/SilentStructFactory.hh>

// AUTO-REMOVED #include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <protocols/topology_broker/TopologyBroker.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/checkpoint/CheckPointer.hh>

#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
//#include <core/scoring/MembraneTopology.fwd.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
//#include <core/pack/task/PackerTask.fwd.hh>

// AUTO-REMOVED #include <protocols/abinitio/Protocol.hh>
// AUTO-REMOVED #include <protocols/moves/TrialMover.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <basic/prof.hh>
#include <utility/exit.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
#include <vector>

#include <protocols/moves/MonteCarlo.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

/// Move these forward declarations to FragmentSampler.fwd.hh
class FragmentSampler;
typedef utility::pointer::owning_ptr< FragmentSampler > FragmentSamplerOP;

//@ brief The Classic Abinitio protocol from rosetta++
/*!
@ detail
general usage:
FragmentSampler  abinitio;
abinitio.init( pose );
...
while(nstruct) {
	 abinitio.apply( pose );
}

call FragmentSampler::register_options() before core::init::init to add relevant options to the applications help

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


class FragmentSampler :  public moves::Mover {    //when changing to new jobdist ====> public moves::Mover {
	typedef moves::Mover Parent;
	typedef moves::Mover BaseClass; //happens to be same as Parent

public:
	virtual ~FragmentSampler();
	//@brief register cmd-line options in option system ( call before core::init::init )
	static void register_options();

	///@brief This constructor does not work -- Fix it before using it.
	// constructor: supply mover classes for Fragment Moves
	FragmentSampler( topology_broker::TopologyBrokerOP broker );

	//@brief FragmentSampler has virtual functions... use this to obtain a new instance
	virtual
	moves::MoverOP clone() const;

	//@brief run protocol on pose
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	//@brief override cycle setting ( init() -> sets it according to cmd-line options )
	virtual void set_cycles( core::Real increase_cycles = 1.0 );

	//@brief set weight - effective for all scoring functions  stage == -1 --> set weight for all stages
	// mod -1 (ignored) ,  in stage3 mod = 1 --> 3a , mod = 2 --> 3b
	virtual void set_score_weight( core::scoring::ScoreType, core::Real setting, StageID stage = ALL_STAGES );

	virtual checkpoint::CheckPointer &get_checkpoints() { return checkpoints_; }

	//	void output_debug_structure( core::pose::Pose&, std::string file_tag ); //make part of Mover class

	void topology_broker( topology_broker::TopologyBrokerOP set );

	//@brief currently used score function ( depends on stage ) -- publich
	core::scoring::ScoreFunction const& current_scorefxn() const;

//	///@brief get membrane topology information
//	core::scoring::MembraneTopologyCOP get_membrane_topology_from_pose(core::pose::Pose const& pose);

protected:
	topology_broker::TopologyBroker const& topology_broker();

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

	//@brief called by constructor ---  calls all set_default_XXX methods
	void set_defaults();

		//@brief set current scorefunction
	void current_scorefxn( core::scoring::ScoreFunction const& scorefxn );

	//@brief If appropriate for the stage then recover low mc.
	void recover_low( core::pose::Pose& pose, StageID stage );

	//@brief this is called if the scorefxn is replaced with a new one
	virtual void replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real intra_stage_progress );

	//@brief set individual weight of current scorefunction
	//@details NOTE: does not change the predefined weights, this information is lost each time the score is set
	//  weight will be overwritten with default ->at beginning of stages, each iteration of loop in stage3
	// we could change this behaviour by having a pointer to the current score and making sure that nobody can
	// change the score inside the mc-object (only const accessor )
	void set_current_weight( core::scoring::ScoreType type, core::Real setting );

	//@brief run cycles for different scoring_stages, return number of steps used
	virtual void do_stage1_cycles( core::pose::Pose &pose );
	virtual void do_stage2_cycles( core::pose::Pose &pose );
	virtual void do_stage3_cycles( core::pose::Pose &pose );
	virtual void do_stage4_cycles( core::pose::Pose &pose );

	//@brief returns true if pose is < 3.0 A rms to last pose sent to this function
	//	bool convergence_check( core::pose::Pose const & pose );

	virtual moves::MoverOP
	mover( core::pose::Pose const& pose, StageID stage_id, core::scoring::ScoreFunction const& scorefxn, core::Real progress = 1.0 );
// 	//@brief returns the Mover that is applied inside the stage1 loop
// 	virtual moves::MoverOP
// 	stage1_mover( core::pose::Pose const& pose );

// 	//@brief returns the Mover that is applied inside the stage3 double loop
// 	virtual moves::MoverOP
// 	stage2_mover( core::pose::Pose const& pose );

// 	//@brief returns the Mover that is applied inside the stage3 double loop
// 	virtual moves::MoverOP
// 	stage3_mover( core::pose::Pose const& pose, core::Size lct1, core::Size lct2 );

// 	//@brief returns the Mover that is applied inside the stage4 loop
// 	virtual moves::MoverOP
// 	stage4_mover( core::pose::Pose const& pose, core::Size kk );


	//@brief called by update_moves() creates the instances of TrialMover with the FragmentMoves
	//virtual void set_trials();

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
	core::Size stage1_cycles() const {
		return stage1_cycles_;
	}

	//@brief Accessor for number of stage2
	inline
	core::Size stage2_cycles() const {
		return stage2_cycles_;
	}

	//@brief Accessor for number of stage3 cycles
	inline
	core::Size stage3_cycles() const {
		return stage3_cycles_;
	}

	//@brief Setter for number of stage4 cycles
	inline
	void set_stage4_cycles( core::Size stage4_cycles_new) {
		stage4_cycles_ = stage4_cycles_new;
	}

	//@brief Accessor for number of stage4 cycles
	inline
	core::Size stage4_cycles() const {
		return stage4_cycles_;
	}

	//@brief query this flag if you have time-intensive stuff and cut it short
	bool bQuickTest() const {
		return bQuickTest_;
	}

	//@brief loophash filter
	bool check_loops(core::pose::Pose& pose);

private:

	// made these set_default methods private since they are not virtual. call set_defaults() to get everything set up
	//@brief read out cmd-line options
	void set_default_options();

	//@brief assigns default score-functions to the 4 stages--> score0 / score1 / (score2/score5) / score3
	void set_default_scores();

	//@brief construct default monto-carlo object
	void set_default_mc( core::scoring::ScoreFunction const& scorefxn	);


	void checkpointed_cycle_block( core::pose::Pose&, StageID, void (FragmentSampler::*cycles) ( core::pose::Pose& ) );

protected:

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarloOP mc_ptr() {
		return mc_;
	}

	//@brief  cycles_numbers  -- modified by set_cycles(), set_default_options()
	Size stage1_cycles_; //score0
	Size stage2_cycles_; //score1
	Size stage3_cycles_; //score25
	Size stage4_cycles_; //score3

private:
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
	bool bQuickTest_;  // this flag might land in base-class

	//@brief a temperature
	core::Real temperature_;

	//@brief a MonteCarlo object -- set_default_mc() , access: mc()
	moves::MonteCarloOP mc_;

	Size total_trials_;

	topology_broker::TopologyBrokerOP topology_broker_;

	checkpoint::CheckPointer checkpoints_;

	utility::vector1< StageID > recover_low_stages_;
	utility::vector1< StageID > skip_stages_;

private:
	static std::string const id2string_[];
	static basic::ProfTag const id2proftag_[];

	std::string const& id2string( StageID id ) {
		assert( id < LAST_STAGE );
		return id2string_[ id ];
	}

	basic::ProfTag id2proftag( StageID id ) {
		assert( id < LAST_STAGE );
		return id2proftag_[ id ];
	}
};

} // abinitio
} // protocols

#endif
