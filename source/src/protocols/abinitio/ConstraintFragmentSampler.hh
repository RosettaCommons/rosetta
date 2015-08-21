// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/ConstraintFragmentSampler.hh
/// @brief header file for ConstraintFragmentSampler protocol
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka

#ifndef INCLUDED_protocols_abinitio_ConstraintFragmentSampler_hh
#define INCLUDED_protocols_abinitio_ConstraintFragmentSampler_hh

// Unit Headers
#include <protocols/abinitio/FragmentSampler.hh>

// Project Headers
#include <protocols/topology_broker/TopologyBroker.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>
#include <vector>

#include <protocols/constraints_additional/MaxSeqSepConstraintSet.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace abinitio {

/// Move these forward declarations to ConstraintFragmentSampler.fwd.hh
class ConstraintFragmentSampler;
typedef utility::pointer::shared_ptr< ConstraintFragmentSampler > ConstraintFragmentSamplerOP;

//@ brief The Classic Abinitio protocol from rosetta++
/*!
@ detail
general usage:
ConstraintFragmentSampler  abinitio;
abinitio.init( pose );
...
while(nstruct) {
abinitio.apply( pose );
}

call ConstraintFragmentSampler::register_options() before core::init::init to add relevant options to the applications help

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

class ConstraintFragmentSampler : public FragmentSampler {
	typedef FragmentSampler Parent;
	typedef moves::Mover BaseClass; //happens to be same as Parent

public:
	/// @brief This constructor does not work -- Fix it before using it.
	// constructor: supply mover classes for Fragment Moves
	ConstraintFragmentSampler( topology_broker::TopologyBrokerOP broker );

	//@brief ConstraintFragmentSampler has virtual functions... use this to obtain a new instance
	virtual
	moves::MoverOP clone() const;

	//@brief run protocol on pose
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	//@brief register cmd-line options in option system ( call before core::init::init )
	static void register_options();

protected:
	/// @brief
	void set_constraint_weight( core::Real setting ) {
		constraint_weight_ = setting;
		set_score_weight( core::scoring::atom_pair_constraint, constraint_weight_ );
	};

	//@brief read out cmd-line options
	void set_defaults();

	//@brief run cycles for different scoring_stages, return number of steps used
	virtual void do_stage1_cycles( core::pose::Pose &pose );

	// anything you want to have done before the stages ?
	//@brief prepare_stageX is called before do_stageX_cycles... overload to change status/scoring/conformation....
	virtual void prepare_stage1( core::pose::Pose &pose );
	virtual void prepare_stage2( core::pose::Pose &pose );

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

	void set_show_viol_level( core::Size setting ) {
		show_viol_level_ = setting;
	}

	virtual void set_max_seq_sep( core::pose::Pose& pose, Size setting );

	core::Real max_seq_sep_fudge() const {
		return max_seq_sep_fudge_;
	}

	void max_seq_sep_fudge( core::Real setting ) {
		max_seq_sep_fudge_ = setting;
	}

	constraints_additional::MaxSeqSepConstraintSet const& constraints() {
		return *constraints_;
	}

	Size total_res( core::pose::Pose const& pose ) const;

	void set_seq_sep_stage1 ( core::Real setting ) {
		seq_sep_stage1_ = setting;
	}

	void set_seq_sep_stage3 ( core::Real setting ) {
		seq_sep_stage3_ = setting;
	}

	virtual void replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real intra_stage_progress );

private:
	core::Real
	evaluate_constraint_energy( core::pose::Pose& pose, core::scoring::ScoreFunction const& ) const;

	constraints_additional::MaxSeqSepConstraintSetOP constraints_;
	core::Real constraint_weight_;

	bool bMinTrial_;

	core::Real max_seq_sep_fudge_;
	core::Real seq_sep_stage1_;
	core::Real seq_sep_stage3_;
	core::Real seq_sep_stage4_;

	//@brief skip cycles in stage1 if nothing is violted ( threshold = ? )
	bool bSkipOnNoViolation_;

	//@brief just for screen output: how verbose should it be
	Size show_viol_level_;

	//@brief usually we do a recover_low before we increase the number of active constraints
	bool bNoRecoverLowAtSwitch_;

	//@brief chainbreak weight is ramped up during stage3 and stage4
	bool bRampChainbreaks_;

	//@brief chainbreak weight is ramped up during stage3 and stage4
	bool bRampCoordConstraints_;

	//@brief overlap chainbreak will be ramped in in stage4
	bool bOverlapChainbreaks_;
};

} // abinitio
} // protocols

#endif
