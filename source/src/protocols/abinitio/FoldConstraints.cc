// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FoldConstraints.cc
/// @brief ab-initio fragment assembly protocol for proteins under the influence of contraints (e.g., NOE)
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/abinitio/FoldConstraints.hh>

// Package headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/numeric.functions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>


//// C++ headers
#include <cstdlib>
#include <string>


// option key includes


#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio.foldconstraints", basic::t_info );

using core::scoring::constraints::ConstraintSet;
using core::scoring::constraints::ConstraintSetOP;
using core::kinematics::ShortestPathInFoldTree;

/// Why is this not being registered  ? Probably not actually needed

void protocols::abinitio::FoldConstraints::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Parent::register_options();
	option.add_relevant( constraints::cst_weight );
	option.add_relevant( fold_cst::no_minimize );
	option.add_relevant( fold_cst::force_minimize );
	option.add_relevant( fold_cst::seq_sep_stages );
	option.add_relevant( fold_cst::reramp_cst_cycles );
	option.add_relevant( fold_cst::reramp_start_cstweight );
	option.add_relevant( fold_cst::reramp_iterations );
	option.add_relevant( fold_cst::skip_on_noviolation_in_stage1 );
	option.add_relevant( fold_cst::stage1_ramp_cst_cycle_factor );
	option.add_relevant( fold_cst::stage2_constraint_threshold );
	option.add_relevant( fold_cst::ignore_sequence_seperation );
	option.add_relevant( fold_cst::no_recover_low_at_constraint_switch );
}

namespace protocols {
namespace abinitio {

using namespace core;

/// @brief c'stor from Movers
FoldConstraints::FoldConstraints(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small,
	int dummy /* otherwise the two constructors are ambigous */
) : ClassicAbinitio( brute_move_small, brute_move_large, smooth_move_small, dummy ),
	constraint_weight_( 1.0 ), run_( 0 )
{
	BaseClass::type( "FoldConstraints" );
	set_default_options();
}

/// @brief c'stor from FragSets --- ClassicFragmentMover and SmoothFragmentMover will be created
FoldConstraints::FoldConstraints(
	core::fragment::FragSetCOP fragset3mer,
	core::fragment::FragSetCOP fragset9mer,
	core::kinematics::MoveMapCOP movemap
) : ClassicAbinitio( fragset3mer, fragset9mer, movemap ),
	constraint_weight_( 1.0 ), run_( 0 )
{
	BaseClass::type( "FoldConstraints" );
	set_default_options();
}

/// @details SHALLOW copy.
FoldConstraints::FoldConstraints( FoldConstraints const & src ) :
	//utility::pointer::ReferenceCount(),
	Parent( src )
{
	min_move_ = src.min_move_;
	constraints_ = src.constraints_;
	constraint_weight_ = src.constraint_weight_;
	bMinTrial_ = src.bMinTrial_;
	bIgnoreSequenceSeparation_ = src.bIgnoreSequenceSeparation_;
	run_ = src.run_;
	max_seq_sep_fudge_ = src.max_seq_sep_fudge_;
	seq_sep_stage1_ = src.seq_sep_stage1_;
	seq_sep_stage3_ = src.seq_sep_stage3_;
	seq_sep_stage4_ = src.seq_sep_stage4_;
	start_ramp_cstweight_ = src.start_ramp_cstweight_;
	ramp_cst_cycles_ = src.ramp_cst_cycles_;
	ramp_iterations_ = src.ramp_iterations_;
	bSkipOnNoViolation_ = src.bSkipOnNoViolation_;
	show_viol_level_ = src.show_viol_level_;
	constraint_threshold_ = src.constraint_threshold_;
}

FoldConstraints::~FoldConstraints() = default;


moves::MoverOP
FoldConstraints::clone() const {
	return moves::MoverOP( new FoldConstraints(*this) );
}


void FoldConstraints::set_default_options() {
	Parent::set_default_options();
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	set_constraint_weight( option[ OptionKeys::constraints::cst_weight ] );
	bMinTrial_ = !bQuickTest() && !option[ OptionKeys::fold_cst::no_minimize ];
	bIgnoreSequenceSeparation_ = option[ OptionKeys::fold_cst::ignore_sequence_seperation ];
	max_seq_sep_fudge_ = 1.0; //models roseta++ behaviour.
	if ( option[ OptionKeys::fold_cst::seq_sep_stages ].user() ) {
		if ( option[ OptionKeys::fold_cst::seq_sep_stages ]().size() != 3 ) {
			utility_exit_with_message("option seq_sep_stages requires exact 3 values!!!");
		}
		seq_sep_stage1_ = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 1 ]; //default 15
		seq_sep_stage3_ = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 2 ]; //default 15
		seq_sep_stage4_ = option[ OptionKeys::fold_cst::seq_sep_stages ]()[ 3 ]; //default 15
	} else {
		seq_sep_stage1_ = 0.15;
		seq_sep_stage3_ = 0.5;
		seq_sep_stage4_ = 1.0;
	}

	start_ramp_cstweight_ = option[ fold_cst::reramp_start_cstweight ];
	ramp_cst_cycles_ =  option[ fold_cst::reramp_cst_cycles ];
	bSkipOnNoViolation_ = option[ fold_cst::skip_on_noviolation_in_stage1 ];
	constraint_threshold_ = option[ fold_cst::stage2_constraint_threshold ];
}

//otherwise stage2 cycles remain as in the classic protocol
Size FoldConstraints::total_res( core::pose::Pose const& pose ) const {
	return static_cast< Size >( std::min( 1.0*pose.total_residue(), constraints_->largest_possible_sequence_sep( pose ) * max_seq_sep_fudge_ ) );
}

// what's an noe_stage?
Size noe_stage( Size total_res, Real factor ) {
	return static_cast< Size >( std::min( factor, 1.0 ) * total_res );
}

bool
FoldConstraints::prepare_stage1( core::pose::Pose& pose ) {
	Parent::prepare_stage1( pose );
	tr.Debug << "active constraints \n";
	return true;
}

void
FoldConstraints::set_max_seq_sep( core::pose::Pose& pose, Size setting ) {
	using namespace basic::options;
	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		tr.Info << "max_seq_sep: " << setting << std::endl;
		// we really have constraints... lets' get to work
		if ( !bIgnoreSequenceSeparation() ) {
			if ( !option[ OptionKeys::fold_cst::no_recover_low_at_constraint_switch ]() ) {
				mc().recover_low( pose ); //is this really what one should do at this point ?
			}
			constraints_->set_max_seq_sep( setting );
			pose.constraint_set( constraints_ );
			if ( !bMinTrial_ )  mc().reset( pose ); //is also done by min_trial... don't do it twice
		}
	} else if ( !option[ OptionKeys::fold_cst::force_minimize ] ) return; //jump out if no constraints and no force_minimize
	if ( bMinTrial_ ) {
		min_trial( pose );
	}
}

bool
FoldConstraints::do_stage1_cycles( pose::Pose& pose ) {
	set_max_seq_sep( pose, std::min(3, int(noe_stage( total_res( pose ), seq_sep_stage1_) ) ) );
	//set_max_seq_sep( pose, 3 ); //yip this is the classic value
	core::Real old_constraint_score = evaluate_constraint_energy ( pose, mc().score_function() );
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	moves::MoverOP trial( stage1_mover( pose, trial_large() ) );
	core::Real const cycle_factor( option[ fold_cst::stage1_ramp_cst_cycle_factor ] );
	Size const cycles ( static_cast< Size > ( cycle_factor * stage1_cycles() ) );
	int total_cycles = 0;
	if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
	//first run a normal set of fragment insertions until extended chain is lost
	total_cycles += Parent::do_stage1_cycles( pose );
	if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );

	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		// Now ramp up the seq_sep of the constraints... still on score0
		for ( Size jk = 3; jk <= noe_stage( total_res( pose ), seq_sep_stage1_); jk += 2 ) {
			//  mc().recover_low( pose ); superfluous -- done in set_max_seq_sep if constraints are actually present
			set_max_seq_sep( pose, jk);
			if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
			if ( old_constraint_score == evaluate_constraint_energy ( pose, mc().score_function() ) ) continue;
			for ( Size j = 1; j <= cycles; ++j, ++total_cycles ) {
				//   if ( evaluate_constraint_energy( pose, mc().score_function() ) < 10.0 ) break; this is unlikely to be triggered for cnc and always triggered for james-cst
				if ( numeric::mod( j, (Size)10)==0 && bSkipOnNoViolation_ && pose.constraint_set()->show_violations( tr, pose, 0 ) == 0 ) break;
				trial->apply( pose );
			}
			old_constraint_score = evaluate_constraint_energy ( pose, mc().score_function() );
		}
	}
	return true;
}

bool
FoldConstraints::prepare_stage2( core::pose::Pose& pose ) {
	set_max_seq_sep( pose, noe_stage( total_res( pose ), seq_sep_stage1_ ) );
	Parent::prepare_stage2( pose );
	if ( tr.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
	// if ( bMinTrial_ ) min_trial( pose ); done in set_max_seq_sep now
	return true;
}

bool
FoldConstraints::do_stage2_cycles( pose::Pose& pose ) {
	set_max_seq_sep( pose, noe_stage( total_res( pose ), seq_sep_stage1_ ) );

	// first do normal stage2 sampling
	bool success = Parent::do_stage2_cycles( pose );

	// if ramp_cycles drop atom_pair_constraint weight to low value and ramp up
	if ( success && ramp_cst_cycles_ > 500 ) {
		moves::MoverOP trial( stage2_mover( pose, trial_large() ) );
		for ( Size loop = 1; loop <= ramp_iterations_; loop++ ) {
			mc().recover_low( pose );
			Real const final_weight( current_scorefxn().get_weight( scoring::atom_pair_constraint ) );
			for ( Size j=1; j<=ramp_cst_cycles_; j++ ) {
				if ( numeric::mod( j, (Size) 50 ) == 0 ) {
					Real weight( (start_ramp_cstweight_ + (1.0-start_ramp_cstweight_)*j/ramp_cst_cycles_ )*final_weight );
					set_current_weight( scoring::atom_pair_constraint, weight );
				}
				trial->apply( pose );
			}
		}
	}
	return success;
}

bool
FoldConstraints::prepare_loop_in_stage3( core::pose::Pose &pose, Size loop_iteration, Size total_iterations ) {
	/* stage3 rosetta++
	noe_stage = 15 + (total_residue/2-15)*kk/nloop;
	classical_constraints::BOUNDARY::set_max_seqSep(noe_stage);
	*/
	bool success = Parent::prepare_loop_in_stage3( pose, loop_iteration, total_iterations );

	if ( constraints_ ) {
		//  mc().recover_low( pose ); superfluous -- recover low in set_max_seq_sep if actually necessary
		core::Real const noe_fact ( seq_sep_stage1_ + ( seq_sep_stage3_ - seq_sep_stage1_) * 1.0*loop_iteration/ total_iterations );
		tr.Debug << "current noe_fact is " << noe_fact << std::endl;
		set_max_seq_sep( pose, noe_stage( total_res( pose ), noe_fact) );
		if ( tr.Info.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );

		if ( constraint_threshold_ > 0 ) {
			constraints_->set_max_seq_sep( pose.total_residue() );
			pose.constraint_set( constraints_ );
			scoring::ScoreFunction cst_score;
			cst_score.set_weight( scoring::atom_pair_constraint, 1.0 );
			core::Real total_cst = cst_score( pose );
			if ( total_cst > constraint_threshold_ ) {
				tr.Info << " Constraint threshold violated: " << total_cst << " (limit: " << constraint_threshold_ << " abandon decoy!" << std::endl;
				return false;
			}
		}

	}
	return success;
}

bool
FoldConstraints::prepare_loop_in_stage4( core::pose::Pose &pose, Size loop_iteration, Size total_iterations ) {
	/* stage3 rosetta++
	noe_stage = 15 + (total_residue/2-15)*kk/nloop;
	classical_constraints::BOUNDARY::set_max_seqSep(noe_stage);
	*/
	bool success = Parent::prepare_loop_in_stage4( pose, loop_iteration, total_iterations );
	// if ( bMinTrial_ ) min_trial( pose ); done in set_max_seq_sep now
	if ( constraints_ ) {
		Real const noe_fact ( seq_sep_stage3_ + ( seq_sep_stage4_ - seq_sep_stage3_ ) * 1.0*loop_iteration / total_iterations );
		tr.Debug << "current noe_fact is " << noe_fact << std::endl;
		set_max_seq_sep( pose, noe_stage( total_res( pose ), noe_fact ) );
	}
	if ( tr.Info.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
	return success;
}

bool
FoldConstraints::prepare_stage4( core::pose::Pose& pose ) {
	bool success = Parent::prepare_stage4( pose );
	if ( tr.Info.visible() ) pose.constraint_set()->show_violations( tr, pose, show_viol_level_ );
	return success;
}
/*
bool
FoldConstraints::do_stage2_cycles( core::pose::Pose& pose ) {
}*/

void
FoldConstraints::set_default_scores() {
	using namespace scoring;
	Parent::set_default_scores();
	tr.Debug << "switch constraints on..." << std::endl;
	set_score_weight( atom_pair_constraint, constraint_weight_ );
	// set_score_weight( angle_constraint, constraint_weight_ );
}


void
FoldConstraints::apply( core::pose::Pose & pose ) {
	// setup constraints to allow for MaxSeqSep
	core::scoring::constraints::ConstraintSetOP orig_constraints( nullptr );
	if ( !bIgnoreSequenceSeparation() ) {
		tr.Debug << "introduce MaxSeqSep Filter for constraints \n";
		orig_constraints = pose.constraint_set()->clone();
		constraints_ = constraints_additional::MaxSeqSepConstraintSetOP( new constraints_additional::MaxSeqSepConstraintSet( *orig_constraints, pose.fold_tree() ) );
		constraints_->set_max_seq_sep( pose.total_residue() ); // so it is prepared for stage4.
		pose.constraint_set( constraints_ );
	}

	mc().reset( pose );
	++run_;
	// call apply of original protocol
	Parent::apply( pose );

	// revert constraints to original form
	if ( orig_constraints ) {
		pose.constraint_set( orig_constraints );
	}
}

std::string
FoldConstraints::get_name() const {
	return "FoldConstraints";
}

void
FoldConstraints::setup_default_min_move() {
	tr.Info << "setup basic minmove" << std::endl;
	min_move_ = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::MinMover );
	min_move_->movemap( movemap() );
	min_move_->min_type( "lbfgs_armijo_nonmonotone" );
}

//@detail overload if your extension stores additional moves as member variables
void FoldConstraints::set_movemap ( core::kinematics::MoveMapCOP mm ) {
	Parent::set_movemap( mm );
	if ( min_move_ ) min_move_->movemap( mm );
}


void FoldConstraints::set_min_move( protocols::simple_moves::MinMoverOP mm) {
	min_move_ = mm;
}

void
FoldConstraints::min_trial( core::pose::Pose &pose ) {
	if ( !min_move_ ) setup_default_min_move();
	tr.Debug << "start minimization... " << std::endl;
	current_scorefxn().show( tr.Debug, pose );
	//get currently used score_function...
	min_move_->score_function( current_scorefxn().clone() );
	min_move_->apply( pose );
	mc().reset( pose );
	tr.Debug << "done with minimization! " << std::endl;
	if ( tr.Info.visible() ) {
		tr.Info << "minimized: ";
		pose.constraint_set()->show_violations( tr.Info, pose, show_viol_level_ );
		tr.Info << std::endl;
	}
}

inline core::Real
FoldConstraints::evaluate_constraint_energy(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & scfxn
) const {
	scfxn( pose ); //accumulate energies

	return
		pose.energies().total_energies()[ core::scoring::atom_pair_constraint  ] +
		pose.energies().total_energies()[ core::scoring::coordinate_constraint ] +
		pose.energies().total_energies()[ core::scoring::dihedral_constraint   ] +
		pose.energies().total_energies()[ core::scoring::angle_constraint      ];

}

} //abinitio
} //protocols


