// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka

// Unit Headers
#include <core/types.hh>

// Package Headers
#include <protocols/abinitio/FragmentSampler.hh>

// Project Headers
#include <protocols/topology_broker/TopologyBroker.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/jd2/util.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// C++ headers
#include <cstdlib>
#include <string>

#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <utility/vector1.hh>

static basic::Tracer tr("protocols.abinitio");

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/*!
@detail call this:
FragmentSampler::register_options() before devel::init().
Derived classes that overload this function should also call Parent::register_options()
*/
void protocols::abinitio::FragmentSampler::register_options() {
	Parent::register_options();
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( OptionKeys::abinitio::increase_cycles );
	option.add_relevant( OptionKeys::abinitio::debug );
	option.add_relevant( OptionKeys::abinitio::skip_stages );
	option.add_relevant( OptionKeys::abinitio::skip_convergence_check );
	option.add_relevant( OptionKeys::abinitio::log_frags );
	option.add_relevant( OptionKeys::abinitio::end_bias );
	option.add_relevant( OptionKeys::abinitio::symmetry_residue );
	option.add_relevant( OptionKeys::abinitio::vdw_weight_stage1 );
	option.add_relevant( OptionKeys::abinitio::override_vdw_all_stages );
	option.add_relevant( OptionKeys::abinitio::recover_low_in_stages );
	option.add_relevant( OptionKeys::abinitio::close_chbrk );
}

namespace protocols {
namespace abinitio {

//little helper function
bool contains_stageid( utility::vector1< abinitio::StageID > vec, abinitio::StageID query ) {
	return find( vec.begin(), vec.end(), query) != vec.end();
}

/// @detail  large (stage1/stage2)
/// small(stage2/stage3/stage4)
/// smooth_small ( stage3/stage4)
FragmentSampler::FragmentSampler( topology_broker::TopologyBrokerOP broker )
	: topology_broker_( broker ),
		checkpoints_("FragmentSampler")
{
	BaseClass::type( "FragmentSampler" );

	set_defaults();
}

FragmentSampler::~FragmentSampler() {}

/// @brief FragmentSampler has virtual functions... use this to obtain a new instance
moves::MoverOP
FragmentSampler::clone() const
{
	return new FragmentSampler( *this );
}

void FragmentSampler::checkpointed_cycle_block(
			 core::pose::Pose& pose,
			 StageID stage_id,
			 void (FragmentSampler::*cycles)(core::pose::Pose& )
) {

	// part X ----------------------------------------
	tr.Info <<  "\n===================================================================\n";
	tr.Info <<  "   " << id2string( stage_id ) << "                                                        \n";
	tr.Info <<  "--------------------------------------------------------------------\n";
	tr.Info << "FragmentSampler: " << id2string( stage_id ) << std::endl;
	PROF_START( id2proftag( stage_id ) );
	clock_t starttime = clock();

	if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), id2string( stage_id ), false /* fullatom*/, true /*fold tree */ )) {
		scoring::constraints::ConstraintSetOP orig_constraints( pose.constraint_set()->clone() );

		//run the fragment cycles
		(this->*cycles)( pose ); //calls do_stageX_cycles()

		recover_low( pose, stage_id );
		if ( tr.Info.visible() ) current_scorefxn().show( tr.Info, pose );
		tr.Info << std::endl;
		mc().show_counters();
		total_trials_+=mc().total_trials();
		mc().reset_counters();

		pose.constraint_set( orig_constraints ); // restore constraints - this is critical for checkpointing to work
		get_checkpoints().checkpoint( pose, get_current_tag(), id2string( stage_id ), true /*fold tree */ );
	} //recover checkpoint
	get_checkpoints().debug( get_current_tag(), id2string( stage_id ), current_scorefxn()( pose ) );

	clock_t endtime = clock();
	PROF_STOP( id2proftag( stage_id ) );
	if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();

	if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
		tr.Info << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
	  jd2::output_intermediate_pose( pose, id2string( stage_id ) );
	}
}

void FragmentSampler::apply( pose::Pose & pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;

	tr.Info << "Fragment Sampler: " << get_current_tag() << std::endl;

	runtime_assert( topology_broker_ ); // really this protocol doesn't make much sense without it
	mc().clear_poses(); // these two statements were only necessary after march 18 2009... something ALF did recently ?
	mc().reset( pose );
	//	current_scorefxn()( pose );
	if ( option[ OptionKeys::run::dry_run ]() ) {
		replace_scorefxn( pose, STAGE_4, 1.0 );
		return;
	}
    
	total_trials_ = 0;
	current_scorefxn()(pose);
  if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
	jd2::output_intermediate_pose( pose, "stage0" );
	}
    
	try {
		if ( !contains_stageid( skip_stages_, STAGE_1 ) ) {
			prepare_stage1( pose );
			checkpointed_cycle_block( pose, STAGE_1, &FragmentSampler::do_stage1_cycles );
		} //skipStage1

		if ( !contains_stageid( skip_stages_, STAGE_2 ) ) {
			prepare_stage2( pose );
			checkpointed_cycle_block( pose, STAGE_2, &FragmentSampler::do_stage2_cycles );
		}

		if ( !contains_stageid( skip_stages_, STAGE_3 ) ) {
			prepare_stage3( pose );
			checkpointed_cycle_block( pose, STAGE_3b, &FragmentSampler::do_stage3_cycles );
		}

		if ( !contains_stageid( skip_stages_, STAGE_4 ) ) {
			prepare_stage4( pose );
			checkpointed_cycle_block( pose, STAGE_4, &FragmentSampler::do_stage4_cycles );
		}

		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Finished Abinitio                                                 \n";
		tr.Info <<  std::endl;
		topology_broker().apply_filter( pose, END_ABINITIO, 1 );
	} catch ( moves::EXCN_Converged& excn ) {
		//		Size last_stage( STAGE_4 );
		//		while( contains_stageid( skip_stages_, last_stage ) ) --last_stage;
		mc().recover_low( pose );
		get_checkpoints().flush_checkpoints();
	};

	replace_scorefxn( pose, STAGE_4, 1.0 ); ///here so we activate all chainbreaks and constraints
	return;
}// FragmentSampler::apply( pose::Pose & pose )


std::string
FragmentSampler::get_name() const {
	return "FragmentSampler";
}

//@detail read cmd_line options and set default versions for many protocol members: trials/moves, score-functions, Monte-Carlo
void FragmentSampler::set_defaults() {
	using namespace basic::options;
	temperature_ = 2.0;
	set_default_scores();
	set_default_options();
	set_default_mc( *score_stage4_ );//if we skip all stages.. this is probably the least surprising "last used scorefxn"
}

//@detail sets Monto-Carlo object to default
void FragmentSampler::set_default_mc(
 	scoring::ScoreFunction const & scorefxn
) {
	set_mc( new moves::MonteCarlo( scorefxn, temperature_ ) );
	canonical_sampling::mc_convergence_checks::setup_convergence_checks_from_cmdline( *mc_ );
}

//@detail sets Monto-Carlo object
void FragmentSampler::set_mc( moves::MonteCarloOP mc_in ) {
	mc_ = mc_in;
}

//@detail override cmd-line setting for "increase_cycling"
void FragmentSampler::set_cycles( Real increase_cycles ) {
	stage1_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage2_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage3_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage4_cycles_ = static_cast< int > (4000 * increase_cycles);

	using namespace basic::options;
}

void FragmentSampler::topology_broker( topology_broker::TopologyBrokerOP set ) {
	topology_broker_ = set;
}

topology_broker::TopologyBroker const& FragmentSampler::topology_broker() {
	runtime_assert( topology_broker_ );
	return *topology_broker_;
}

void FragmentSampler::set_default_scores() {
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "creating standard scoring functions" << std::endl;

	if ( option[ OptionKeys::abinitio::stage1_patch ].user() ) {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0", option[ OptionKeys::abinitio::stage1_patch ]() );
	} else {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0" );
	}

	if ( option[ OptionKeys::abinitio::stage2_patch ].user() ) {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1", option[ OptionKeys::abinitio::stage2_patch ]() );
	} else {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1" );
	}

	if ( option[ OptionKeys::abinitio::stage3a_patch ].user() ) {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2", option[ OptionKeys::abinitio::stage3a_patch ]() );
	} else {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2" );
	}

	if ( option[ OptionKeys::abinitio::stage3b_patch ].user() ) {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5", option[ OptionKeys::abinitio::stage3b_patch ]() );
	} else {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5" );
	}

	if ( option[ OptionKeys::abinitio::stage4_patch ].user() ) {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage4_patch ]() );
	} else {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3" );
	}

	if ( option[ OptionKeys::abinitio::override_vdw_all_stages ] ) {
		set_score_weight( scoring::vdw, option[ OptionKeys::abinitio::vdw_weight_stage1 ], ALL_STAGES );
	}
}


/// @brief sets a score weight for all stages of abinitio
void FragmentSampler::set_score_weight( scoring::ScoreType type, Real setting, StageID stage ) {
	tr.Debug << "set score weights for ";
	if ( stage == ALL_STAGES ) tr.Debug << "all stages ";
	else tr.Debug << "stage " << (stage <= STAGE_3a ? stage : ( stage-1 ) ) << ( stage == STAGE_3b ? "b " : " " );
	tr.Debug << scoring::name_from_score_type(type) << " " << setting << std::endl;
	if (score_stage1_  && ( stage == STAGE_1  || stage == ALL_STAGES )) score_stage1_ ->set_weight(type, setting);
	if (score_stage2_  && ( stage == STAGE_2  || stage == ALL_STAGES )) score_stage2_ ->set_weight(type, setting);
	if (score_stage3a_ && ( stage == STAGE_3a || stage == ALL_STAGES )) score_stage3a_->set_weight(type, setting);
	if (score_stage3b_ && ( stage == STAGE_3b || stage == ALL_STAGES )) score_stage3b_->set_weight(type, setting);
	if (score_stage4_  && ( stage == STAGE_4  || stage == ALL_STAGES )) score_stage4_ ->set_weight(type, setting);
}

//@brief currently used score function ( depends on stage )
scoring::ScoreFunction const& FragmentSampler::current_scorefxn() const {
	return mc().score_function();
}

//@brief set current scorefunction
void FragmentSampler::current_scorefxn( scoring::ScoreFunction const& scorefxn ) {
	mc().score_function( scorefxn );
}

//@brief set individual weight of current scorefunction --- does not change the predefined scores: score_stageX_
void FragmentSampler::set_current_weight( core::scoring::ScoreType type, core::Real setting ) {
	scoring::ScoreFunctionCOP scorefxn_const ( mc().score_function() );
	scoring::ScoreFunctionOP scorefxn = scorefxn_const->clone();
	scorefxn->set_weight( type, setting );
	mc().score_function( *scorefxn ); //trigger rescore
}

void FragmentSampler::set_default_options() {
	using namespace basic::options;
	just_smooth_cycles_ = option[ OptionKeys::abinitio::smooth_cycles_only ]; // defaults to false
	bQuickTest_ = basic::options::option[ basic::options::OptionKeys::run::test_cycles ]();

	if ( bQuickTest() ) {
		set_cycles( 0.001 );
	} else {
		set_cycles( option[ OptionKeys::abinitio::increase_cycles ] ); // defaults to factor of 1.0
	}

	apply_large_frags_   = true;  // apply large frags in phase 2!

	// in rosetta++ switched on in fold_abinitio if contig_size < 30 in pose_abinitio never
	short_insert_region_ = false;  // apply small fragments in phase 2!

	skip_stages_.clear();
	if ( option[ OptionKeys::abinitio::skip_stages ].user() ) {
		for ( IntegerVectorOption::const_iterator it = option[ OptionKeys::abinitio::skip_stages ]().begin(),
					eit = option[ OptionKeys::abinitio::skip_stages ]().end(); it!=eit; ++it ) {
			if ( *it == 1 ) skip_stages_.push_back( STAGE_1 );
			else if ( *it == 2 ) skip_stages_.push_back( STAGE_2 );
			else if ( *it == 3 ) skip_stages_.push_back( STAGE_3 );
			else if ( *it == 4 ) skip_stages_.push_back( STAGE_4 );
		}
	}

	if ( option[ OptionKeys::abinitio::recover_low_in_stages ].user() ) {
		for ( IntegerVectorOption::const_iterator it = option[ OptionKeys::abinitio::recover_low_in_stages ]().begin(),
					eit = option[ OptionKeys::abinitio::recover_low_in_stages ]().end(); it!=eit; ++it ) {
			if ( *it == 1 ) recover_low_stages_.push_back( STAGE_1 );
			else if ( *it == 2 ) recover_low_stages_.push_back( STAGE_2 );
			else if ( *it == 3 ) {
				recover_low_stages_.push_back( STAGE_3a );
				recover_low_stages_.push_back( STAGE_3b );
			}
			else if ( *it == 4 ) recover_low_stages_.push_back( STAGE_4 );
		}
	} else {
		recover_low_stages_.clear();
		recover_low_stages_.push_back( STAGE_1 );
		recover_low_stages_.push_back( STAGE_2 );
		recover_low_stages_.push_back( STAGE_3a );
		recover_low_stages_.push_back( STAGE_3b );
		recover_low_stages_.push_back( STAGE_4 );
	}
}

moves::MoverOP
FragmentSampler::mover( pose::Pose const& pose, StageID stage_id, core::scoring::ScoreFunction const& scorefxn, core::Real progress ) {
	topology_broker().apply_filter( pose, stage_id, progress ); //raises exception if filter failed
	return topology_broker().mover( pose, stage_id, scorefxn, progress );
}

void FragmentSampler::do_stage1_cycles( pose::Pose &pose ) {
	moves::RepeatMover( new moves::TrialMover( mover( pose, STAGE_1, current_scorefxn() ), mc_ptr() ), stage1_cycles() ).apply( pose );
	mc().reset( pose ); // make sure that we keep the final structure
  if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
	jd2::output_intermediate_pose( pose, "stage1_cycles" );
	}
}

void FragmentSampler::do_stage2_cycles( pose::Pose &pose ) {
	moves::RepeatMover( new moves::TrialMover( mover( pose, STAGE_2, current_scorefxn() ), mc_ptr() ), stage2_cycles() ).apply( pose );
}

/*! @detail stage3 cycles:
	nloop1 : outer iterations
	nloop2 : inner iterations
	stage3_cycle : trials per inner iteration
	every inner iteration we switch between score_stage3a ( default: score2 ) and score_stage3b ( default: score 5 )

	prepare_loop_in_stage3() is called before the stage3_cycles() of trials are started.

	first outer loop-iteration is done with TrialMover trial_large()
	all following iterations with trial_small()

	start each iteration with the lowest_score_pose. ( mc->recover_low() -- called in prepare_loop_in_stage3() )

*/
void FragmentSampler::do_stage3_cycles( pose::Pose &pose ) {
	using namespace ObjexxFCL;

	// interlaced score2 / score 5 loops
	// nloops1 and nloops2 could become member-variables and thus changeable from the outside
	int nloop1 = 1;
	int nloop2 = 10; //careful: if you change these the number of structures in the structure store changes.. problem with checkpointing
	// individual checkpoints for each stage3 iteration would be a remedy. ...

	if ( short_insert_region_ ) {
		nloop1 = 2;
		nloop2 = 5;
	}
	Size const total_iterations ( nloop1*nloop2 );

	int iteration = 1;
	for ( int lct1 = 1; lct1 <= nloop1; lct1++) {
		for ( int lct2 = 1; lct2 <= nloop2; lct2++, iteration++  ) {
			tr.Debug << "Loop: " << lct1 << "   " << lct2 << std::endl;
			StageID current_stage_id = ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) ? STAGE_3a : STAGE_3b;
			prepare_loop_in_stage3( pose, iteration, nloop1*nloop2 );
			core::Real progress( 1.0* iteration/total_iterations );

			if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2),
			                                       false /*fullatom */, true /*fold tree */ )) {

				tr.Debug << "  Score stage3 loop iteration " << lct1 << " " << lct2 << std::endl;
				moves::TrialMoverOP stage3_trials = new moves::TrialMover( mover( pose, current_stage_id, current_scorefxn(), progress ), mc_ptr() );
					moves::RepeatMover( stage3_trials, stage3_cycles() ).apply( pose );

				recover_low( pose, current_stage_id );
				get_checkpoints().checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), true /*fold tree */ );
			}//recover_checkpoint
			get_checkpoints().debug( get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), current_scorefxn()( pose ) );
		} // loop 2
	} // loop 1
}

void FragmentSampler::do_stage4_cycles( pose::Pose &pose ) {
	Size nloop_stage4 = 3;
	for ( Size kk = 1; kk <= nloop_stage4; ++kk ) {
		prepare_loop_in_stage4( pose, kk, nloop_stage4 );

		if (!get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), false /* fullatom */, true /* fold_tree */ )) {
			tr.Debug << "start " << stage4_cycles() << " cycles" << std::endl;
			moves::TrialMoverOP stage4_trials = new moves::TrialMover( mover( pose, STAGE_4, current_scorefxn(), 1.0*kk/nloop_stage4 ), mc_ptr() );
				moves::RepeatMover( stage4_trials, stage4_cycles() ).apply( pose );
			tr.Debug << "finished" << std::endl;

			recover_low( pose, STAGE_4 );
			get_checkpoints().checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), true /*fold tree */ );
		}
		get_checkpoints().debug( get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk),  current_scorefxn()( pose ) );

		//don't store last structure since it will be exactly the same as the final structure delivered back via apply
	}  // loop kk
}

void FragmentSampler::recover_low( core::pose::Pose& pose, StageID stage ){
	if ( contains_stageid( recover_low_stages_, stage ) ) {
		if ( tr.Trace.visible() ) current_scorefxn().show( tr.Trace, pose );
		tr.Trace << "================ RECOVER LOW ==============" << std::endl;
		mc_->recover_low( pose );
		if ( tr.Trace.visible() ) current_scorefxn().show( tr.Trace, pose );
		tr.Trace << std::endl;
	} else {
		tr.Trace << "================ RECOVER LOW SKIPPED ===============" << std::endl;
	}
}

// anything you want to have done before the stages ?
void FragmentSampler::replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real /*intra_stage_progress */ ) {
	// must assume that the current pose is the one to be accepted into the next stage! (this change was necessary for
	// checkpointing to work correctly. --> that means no recover_low at this point! have to call it explicitly in stage3_loops

	//intra_stage_progress = intra_stage_progress;
	if (score_stage1_  && ( stage == STAGE_1 )) current_scorefxn( *score_stage1_ );
	if (score_stage2_  && ( stage == STAGE_2 )) current_scorefxn( *score_stage2_ );
	if (score_stage3a_ && ( stage == STAGE_3a || stage == STAGE_3 )) current_scorefxn( *score_stage3a_ );
	if (score_stage3b_ && ( stage == STAGE_3b)) current_scorefxn( *score_stage3b_ );
	if (score_stage4_  && ( stage == STAGE_4 )) current_scorefxn( *score_stage4_ );
	mc_->set_autotemp( true, temperature_ );
	mc_->set_temperature( temperature_ ); // temperature might have changed due to autotemp..
	mc_->reset( pose );

	current_scorefxn()(pose);
}

// prepare stage1 sampling
void FragmentSampler::prepare_stage1( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_1, 0.5 );
	mc_->set_autotemp( false, temperature_ );
	/// Now handled automatically.  score_stage1_->accumulate_residue_total_energies( pose ); // fix this
}

void FragmentSampler::prepare_stage2( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_2, 0.5 );
}

void FragmentSampler::prepare_stage3( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_3a, 0 );
}

void FragmentSampler::prepare_stage4( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_4, 0 );
}

void FragmentSampler::prepare_loop_in_stage3( core::pose::Pose &pose/*pose*/, Size iteration, Size total ){
	if ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) {
		replace_scorefxn( pose, STAGE_3a, 1.0* iteration/total );
	} else {
		replace_scorefxn( pose, STAGE_3b, 1.0* iteration/total );
	}
}

void FragmentSampler::prepare_loop_in_stage4( core::pose::Pose &pose, Size iteration, Size total ){
	replace_scorefxn( pose, STAGE_4, 1.0* iteration/total );
}

std::string const FragmentSampler::id2string_[] = { "all_stages", "stage1", "stage2", "stage3", "stage3", "stage3", "stage4"};

basic::ProfTag const FragmentSampler::id2proftag_[] = { basic::STAGE1/*ALL_STAGES*/, basic::STAGE1, basic::STAGE2, basic::STAGE3, basic::STAGE3, basic::STAGE3, basic::STAGE4 };

} //abinitio
} //protocols
