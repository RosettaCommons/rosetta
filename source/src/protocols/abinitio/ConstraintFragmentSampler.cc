// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ConstraintFragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/abinitio/ConstraintFragmentSampler.hh>
#include <protocols/topology_broker/TopologyBroker.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/moves/TrialMover.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// C++ headers
#include <cstdlib>
#include <string>

// Auto Headers
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>
// AUTO-REMOVED #include <utility/io/mpistream.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.abinitio" );

using core::Real;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::methods::EnergyMethodOptions;
using std::string;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/*!
@detail call this:
ConstraintFragmentSampler::register_options() before devel::init().
Derived classes that overload this function should also call Parent::register_options()
*/

Size seq_sep_stage( core::Size total_res, Real factor ) {
	return static_cast< Size >( std::min( factor, 1.0 ) * total_res );
}

void protocols::abinitio::ConstraintFragmentSampler::register_options() {
	Parent::register_options();
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( constraints::cst_weight );
	option.add_relevant( fold_cst::no_minimize );
	option.add_relevant( fold_cst::force_minimize );
	option.add_relevant( fold_cst::seq_sep_stages );
	option.add_relevant( fold_cst::skip_on_noviolation_in_stage1 );
	option.add_relevant( fold_cst::stage1_ramp_cst_cycle_factor );
	option.add_relevant( fold_cst::ignore_sequence_seperation );
	option.add_relevant( fold_cst::no_recover_low_at_constraint_switch );
	option.add_relevant( fold_cst::ramp_coord_cst );
}

namespace protocols {
namespace abinitio {

/// @detail  large (stage1/stage2)
/// small(stage2/stage3/stage4)
/// smooth_small ( stage3/stage4)
ConstraintFragmentSampler::ConstraintFragmentSampler(topology_broker::TopologyBrokerOP broker)
	: Parent(broker) {
	BaseClass::type( "ConstraintFragmentSampler" );
	set_defaults();
}

/// @brief ConstraintFragmentSampler has virtual functions... use this to obtain a new instance
moves::MoverOP ConstraintFragmentSampler::clone() const {
	return new ConstraintFragmentSampler( *this );
}

//@detail read cmd_line options and set default versions for many protocol members: trials/moves, score-functions, Monte-Carlo
void ConstraintFragmentSampler::set_defaults() {
	Parent::set_defaults();
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ OptionKeys::constraints::cst_weight ].user() ){
		set_constraint_weight( option[ OptionKeys::constraints::cst_weight ] );

	}
	show_viol_level_ = option[OptionKeys::constraints::viol_level];
	max_seq_sep_fudge_ = 1.0; //models rosetta++ behaviour.

	if ( option[ OptionKeys::fold_cst::seq_sep_stages ].user() ) {
		if ( option[ OptionKeys::fold_cst::seq_sep_stages ]().size() != 3 ) {
			utility_exit_with_message("option seq_sep_stages requires exact 3 values!!!");
		}
		seq_sep_stage1_ = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 1 ]; //default 15
		seq_sep_stage3_ = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 2 ]; //default 15
		seq_sep_stage4_ = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 3 ]; //default 15
	} else {
		seq_sep_stage1_ = 0.15;
		seq_sep_stage3_ = 1.0;  //have basically all but the furthest chainbreaks active at end of stage3.
		seq_sep_stage4_ = 1.2; //make sure that all chainbreaks are switched on in the end
	}

	bSkipOnNoViolation_ = option[ fold_cst::skip_on_noviolation_in_stage1 ];

	bNoRecoverLowAtSwitch_ = option[ fold_cst::no_recover_low_at_constraint_switch ];

	bRampChainbreaks_ = option[ jumps::ramp_chainbreaks ]; //default is true
	bRampCoordConstraints_ = option[ fold_cst::ramp_coord_cst ]; //default is false
  bOverlapChainbreaks_ = option[ jumps::overlap_chainbreak ];
	Real chainbreak_score_1 = option[ jumps::chainbreak_weight_stage1 ]();
  Real chainbreak_score_2 = option[ jumps::chainbreak_weight_stage2 ]();
  Real chainbreak_score_3 = option[ jumps::chainbreak_weight_stage3 ]();
  Real chainbreak_score_4 = option[ jumps::chainbreak_weight_stage4 ]();
  if ( !bRampChainbreaks_ ) {
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_1, STAGE_1 );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_2, STAGE_2 );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_3, STAGE_3a );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_3, STAGE_3b );
    set_score_weight( scoring::linear_chainbreak, chainbreak_score_4, STAGE_4 );
  }
}

void ConstraintFragmentSampler::replace_scorefxn(core::pose::Pose& pose,
                                                 StageID stage,
                                                 Real intra_stage_progress ) {
 	Parent::replace_scorefxn( pose, stage, intra_stage_progress );

	if ( bRampChainbreaks_ ) { //should this live in replace_scorefxn ???
    Real setting( 0.0 );

    if ( stage == STAGE_2 ) {
			setting = 0.25 / 3;
		} else if ( stage == STAGE_3a ) {
			setting = 2.5 * intra_stage_progress * 1.0/3;
		} else if ( stage == STAGE_3b ) {
			setting = 0.5 * intra_stage_progress * 1.0/3;
		} else if ( stage == STAGE_4 ) {
			setting = (1.5*intra_stage_progress+2.5 ) * ( 1.0/3);
		}
		set_current_weight( scoring::linear_chainbreak, setting * option[ jumps::increase_chainbreak ] ); //default 1.0
		if ( bRampCoordConstraints_ ) {
			set_current_weight( scoring::coordinate_constraint, setting * option[ jumps::increase_chainbreak ] );
		}

		if ( bOverlapChainbreaks_ && stage == STAGE_4 ) {
      set_score_weight( scoring::overlap_chainbreak, intra_stage_progress * option[ jumps::increase_chainbreak ] , STAGE_4 );
			set_current_weight( scoring::overlap_chainbreak, intra_stage_progress * option[ jumps::increase_chainbreak ] );
    }
  }
	Real sep_fact( seq_sep_stage1_ );
	if ( stage == STAGE_1 ) {
		sep_fact = std::min( 3.0/total_res( pose ), seq_sep_stage1_ );
	} else if ( stage == STAGE_2 ) {
		sep_fact = seq_sep_stage1_;
	} else if ( stage == STAGE_3a || stage == STAGE_3b ) {
		sep_fact = seq_sep_stage1_ + ( seq_sep_stage3_	- seq_sep_stage1_) * 1.0*intra_stage_progress;
	} else if ( stage == STAGE_4 ) {
		sep_fact = seq_sep_stage3_ + ( seq_sep_stage4_ - seq_sep_stage3_ ) * 1.0*intra_stage_progress;
	}
	tr.Debug << "replace_scorefxn... sep_fact: " << sep_fact << " STAGE " << stage << std::endl;
	tr.Debug << "total number of residues: " << pose.total_residue() << " max. separation in fold-tree " << total_res(pose) << std::endl;
	set_max_seq_sep( pose, seq_sep_stage( total_res( pose ), sep_fact ) );
}

//otherwise stage2 cycles remain as in the classic protocol
Size ConstraintFragmentSampler::total_res( core::pose::Pose const& pose ) const {
	return static_cast< Size >(
			std::min( 1.0*pose.total_residue(), constraints_->largest_possible_sequence_sep( pose ) * max_seq_sep_fudge_ )
	);
}

void ConstraintFragmentSampler::set_max_seq_sep( pose::Pose& pose, core::Size setting ) {
	using namespace basic::options;
	using core::pose::Pose;
	using core::scoring::ScoreFunction;

	bool const bHaveConstraints( pose.constraint_set()->has_residue_pair_constraints() );
	bool const bHaveChainbreaks( topology_broker().has_chainbreaks_to_close() );

	if (bHaveConstraints)
		runtime_assert( constraints_ != 0 );

	if ( bHaveConstraints || bHaveChainbreaks ) {
		if ( bHaveConstraints) {
			// pass along the updated sequence separation parameter to the constraints
			constraints_->set_max_seq_sep( setting );
		}
		tr.Info << "max_seq_sep: " << setting << std::endl;

		// we really have constraints... let's get to work
		if (!bNoRecoverLowAtSwitch_) {
			// Recover the low energy pose, otherwise we will lose it at the next reset()
			if (tr.Trace.visible())
				current_scorefxn().show(tr.Trace, pose);

			tr.Trace << "================ RECOVER LOW (max_seq_sep) ==============" << std::endl;
			mc().recover_low(pose);

			if (bHaveConstraints)
				pose.constraint_set(constraints_);

			if (bHaveChainbreaks) {
				// Take note of our existing energy method options
				const EnergyMethodOptions& current_options =
					current_scorefxn().energy_method_options();

				// Update the maximum sequence separation setting, retaining the values
				// of all other energy method options
				EnergyMethodOptions updated_options(current_options);
				updated_options.cst_max_seq_sep(setting);

				// Replace the score function
                core::scoring::ScoreFunctionOP new_scorefxn = current_scorefxn().clone();
				new_scorefxn->set_energy_method_options(updated_options);
				current_scorefxn(*new_scorefxn);
			}

			mc().reset(pose);
			if (tr.Trace.visible())
				current_scorefxn().show(tr.Trace, pose);
			tr.Trace << std::endl;
		} else {
			Pose low_pose;
			mc().recover_low(low_pose);

			if (bHaveConstraints)
				low_pose.constraint_set(constraints_);

			if (bHaveChainbreaks) {
        // Take note of our existing energy method options
        const EnergyMethodOptions& current_options =
          current_scorefxn().energy_method_options();

        // Update the maximum sequence separation setting, retaining the values
        // of all other energy method options
        EnergyMethodOptions updated_options(current_options);
        updated_options.cst_max_seq_sep(setting);

        // Replace the score function
        ScoreFunctionOP new_scorefxn(current_scorefxn().clone());
        new_scorefxn->set_energy_method_options(updated_options);
        current_scorefxn(*new_scorefxn);
			}

			mc().reset(low_pose);

			if (bHaveConstraints)
				pose.constraint_set(constraints_);

			if (bHaveChainbreaks) {
        // Take note of our existing energy method options
        const EnergyMethodOptions& current_options =
          current_scorefxn().energy_method_options();

        // Update the maximum sequence separation setting, retaining the values
        // of all other energy method options
        EnergyMethodOptions updated_options(current_options);
        updated_options.cst_max_seq_sep(setting);

        // Replace the score function
        ScoreFunctionOP new_scorefxn(current_scorefxn().clone());
        new_scorefxn->set_energy_method_options(updated_options);
        current_scorefxn(*new_scorefxn);
			}

			mc().boltzmann(pose, "add_constraints");
		}
	}
}

void ConstraintFragmentSampler::do_stage1_cycles( pose::Pose& pose ) {
	Real old_constraint_score = evaluate_constraint_energy ( pose, mc().score_function() );

	moves::MoverOP trial( new moves::TrialMover( mover( pose, STAGE_1, current_scorefxn()  ), mc_ptr() ) );

	if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );

	//first run a normal set of fragment insertions until extended chain is lost
	Parent::do_stage1_cycles( pose );

	if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );

	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		// Now ramp up the seq_sep of the constraints... still on score0
		for ( Size jk = 3; jk <= seq_sep_stage( total_res( pose ), seq_sep_stage1_); jk += 2 ) {
			set_max_seq_sep( pose, jk);
			if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
			if ( std::abs( old_constraint_score - evaluate_constraint_energy ( pose, mc().score_function() ) ) < 0.01 ) continue;
			for ( Size j = 1; j <= stage1_cycles(); ++j ) {
				if ( numeric::mod( j, (Size) 10 ) == 0 && bSkipOnNoViolation_ && pose.constraint_set()->show_violations( tr, pose, 0 ) == 0 ) break;
				trial->apply( pose );
			}
			old_constraint_score = evaluate_constraint_energy ( pose, mc().score_function() );
		}
	}
}

void
ConstraintFragmentSampler::prepare_stage1( core::pose::Pose &pose ) {
	//we set score term in here: we minimize in this one... need to set all score terms before that
	Parent::prepare_stage1( pose );
}

void
ConstraintFragmentSampler::prepare_stage2( core::pose::Pose& pose ) {
	//	set_max_seq_sep( pose, seq_sep_stage( total_res( pose ), seq_sep_stage1_ ) );
	Parent::prepare_stage2( pose );
	if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
}

void
ConstraintFragmentSampler::prepare_loop_in_stage3( core::pose::Pose &pose, Size loop_iteration, Size total_iterations ) {
	Parent::prepare_loop_in_stage3( pose, loop_iteration, total_iterations );
	if ( tr.Info.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
}

void
ConstraintFragmentSampler::prepare_loop_in_stage4( core::pose::Pose &pose, Size loop_iteration, Size total_iterations ) {
	Parent::prepare_loop_in_stage4( pose, loop_iteration, total_iterations );
	if ( tr.Info.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
}

void ConstraintFragmentSampler::apply(core::pose::Pose& pose) {
	using core::scoring::constraints::ConstraintSetOP;
	tr.Info << "ConstraintFragment Sampler: " << get_current_tag() << std::endl;

  // take note of the current constraints, as we are responsible for restoring them
  ConstraintSetOP orig_constraints(NULL);
	tr.Debug << "introduce MaxSeqSep Filter for constraints" << std::endl;
	orig_constraints = pose.constraint_set()->clone();

  // initialize a MaxSeqSepConstraintSet with the current set of constraints
	constraints_ = new constraints_additional::MaxSeqSepConstraintSet(*orig_constraints, pose.fold_tree());
	constraints_->set_max_seq_sep(pose.total_residue()); // so it is prepared for stage4.

  // replace <pose>'s ConstraintSet with our newly initialized MaxSeqSepConstraintSet
	pose.constraint_set(constraints_);

	mc().clear_poses();
	mc().reset(pose);

  // *** Parent::apply() ***
	Parent::apply(pose);

  // restore the original set of constraints
	if (orig_constraints)
		pose.constraint_set(orig_constraints);
}

string ConstraintFragmentSampler::get_name() const {
	return "ConstraintFragmentSampler";
}

inline Real ConstraintFragmentSampler::evaluate_constraint_energy(core::pose::Pose& pose,
	                                                                ScoreFunction const & scfxn) const {
	scfxn( pose ); //accumulate energies
	return
		pose.energies().total_energies()[ core::scoring::atom_pair_constraint  ] +
		pose.energies().total_energies()[ core::scoring::coordinate_constraint ] +
		pose.energies().total_energies()[ core::scoring::dihedral_constraint   ] +
		pose.energies().total_energies()[ core::scoring::angle_constraint      ];
}

} //abinitio
} //protocols
