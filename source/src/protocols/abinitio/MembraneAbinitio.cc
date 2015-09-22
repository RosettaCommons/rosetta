// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MembraneAbinitio.cc
/// @brief ab-initio fragment assembly protocol for membrane proteins
/// @details
///   Contains currently: Membrane Abinitio
///
///
/// @author Bjorn Wallner (copied some time ago from ClassicAbinitio.cc of Oliver Lange )


// Unit Headers
#include <protocols/abinitio/MembraneAbinitio.hh>
#include <protocols/simple_moves/GunnCost.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembranePotential.fwd.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <basic/datacache/BasicDataCache.hh>

// Option headers
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/WhileMover.hh>

#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif


// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#ifdef WIN32
#include <ctime>
#endif

//Auto Headers
#include <core/fragment/FragSet.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


static THREAD_LOCAL basic::Tracer tr( "protocols.membrane.abinitio", basic::t_info );


using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/*!
@detail call this:
MembraneAbinitio::register_options() before devel::init().
Derived classes that overload this function should also call Parent::register_options()
*/

OPT_1GRP_KEY( Boolean, abinitio, debug_structures_2 )
//OPT_1GRP_KEY( Boolean, abinitio, membrane_print )
OPT_1GRP_KEY( Boolean, abinitio, test_2 )
OPT_1GRP_KEY( File, abinitio, log_frags_2 )
OPT_1GRP_KEY( Boolean, abinitio, only_stage1_2 )
OPT_1GRP_KEY( Real, abinitio, end_bias_2 )
OPT_1GRP_KEY( Integer, templates, change_movemap_2 )
OPT_1GRP_KEY( Integer, abinitio, symmetry_residue_2 )

void protocols::abinitio::MembraneAbinitio::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( OptionKeys::abinitio::increase_cycles );
	option.add_relevant( OptionKeys::abinitio::smooth_cycles_only );
	option.add_relevant( OptionKeys::abinitio::debug );
	option.add_relevant( OptionKeys::abinitio::skip_convergence_check );
	NEW_OPT( abinitio::debug_structures_2, "write structures to debug-out files", false );
	//NEW_OPT( abinitio::membrane_print, "be noisy and very informative...", false );
	NEW_OPT( abinitio::test_2, "if set the protocol will run quickly to test your setup (files, cmd-line, etc. )", false );
	NEW_OPT( abinitio::log_frags_2, "fragment insertions (each trial) will be logged to file", "" );
	NEW_OPT( abinitio::only_stage1_2, "useful for benchmarks sets cycle of all higher stages to 0", false );
	NEW_OPT( abinitio::end_bias_2, "set the endbias for Fragment moves", 30.0 );
	NEW_OPT( templates::change_movemap_2, "stage in which movemap is switched to allow all bb-residues to move, valid stages: 3..4 (HACK)", 3);
	NEW_OPT( abinitio::symmetry_residue_2,"hacky symmetry mode for dimers, fragments are inserted at i and i + SR - 1", -1 );
}


namespace protocols {
namespace abinitio {

/// @detail  large (stage1/stage2)
/// small(stage2/stage3/stage4)
/// smooth_small ( stage3/stage4)
MembraneAbinitio::MembraneAbinitio(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_small_top25,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small,
	int  /*dummy otherwise the two constructors are ambiguous */
) :
	brute_move_small_( brute_move_small ),
	brute_move_small_top25_( brute_move_small_top25 ),
	brute_move_large_( brute_move_large ),
	smooth_move_small_( smooth_move_small ),
	checkpoint_( "MembraneAbinitio" )
	// output_tag_( "debug ")
{
	BaseClass::type( "MembraneAbintio" );

	// std::cerr << "MembraneAbinitio::constructor has stubbed out...(fatal) see code file";
	// assert( 0 ); //---> needs more implementation to use this constructor: e.g. read out movemap from FragmentMover...
	movemap_ = brute_move_large->movemap();
	//  set_defaults( pose ); in constructor virtual functions are not called
}

MembraneAbinitio::MembraneAbinitio(
	core::fragment::FragSetCOP fragset_small,
	core::fragment::FragSetCOP fragset_small_top25,
	core::fragment::FragSetCOP fragset_large,
	core::kinematics::MoveMapCOP movemap
)  :
	movemap_( movemap ),
	checkpoint_( "MembraneAbinitio" )
{
	BaseClass::type( "MembraneAbinitio" );
	using namespace basic::options;
	simple_moves::ClassicFragmentMoverOP bms, bms25, bml, sms;
	/* if ( option[ OptionKeys::abinitio::log_frags_2 ].user() ) {
	if ( !option[ OptionKeys::abinitio::debug ] ) utility_exit_with_message( "apply option abinitio::log_frags always together with abinitio::debug!!!");
	bms  = new LoggedFragmentMover( fragset_small, movemap );
	bml  = new LoggedFragmentMover( fragset_large, movemap );
	sms  = new SmoothFragmentMover( fragset_small, movemap, new GunnCost );
	} else if ( option[ OptionKeys::abinitio::symmetry_residue_2 ].user() ) {
	Size const sr (  option[ OptionKeys::abinitio::symmetry_residue_2 ] );
	bms = new SymmetricFragmentMover( fragset_small, movemap, sr );
	bml = new SymmetricFragmentMover( fragset_large, movemap, sr );
	sms = new SmoothSymmetricFragmentMover( fragset_small, movemap, new GunnCost, sr );
	} else { */

	using namespace protocols::simple_moves;
	bms25 = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_small_top25, movemap ) );
	bms = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_small, movemap ) );
	bml = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_large, movemap ) );
	sms = simple_moves::ClassicFragmentMoverOP( new SmoothFragmentMover ( fragset_small, movemap, FragmentCostOP( new GunnCost ) ) );
	// }
	/*
	bms->set_end_bias( option[ OptionKeys::abinitio::end_bias_2 ] ); //default is 30.0
	bml->set_end_bias( option[ OptionKeys::abinitio::end_bias_2 ] );
	bms->set_end_bias( option[ OptionKeys::abinitio::end_bias_2 ] );
	*/

	brute_move_small_ = bms;
	brute_move_small_top25_ = bms25;
	brute_move_large_ = bml;
	smooth_move_small_ = sms;
}

/// @brief setup moves, mc-object, scores
/// @details can't call this from constructor; virtual functions don't operate until construction has completed.

void
MembraneAbinitio::init( core::pose::Pose const& pose ) {
	// Parent::init( pose );
	set_defaults( pose );
	// bInitialized_ = true;
}

/// @brief MembraneAbinitio has virtual functions... use this to obtain a new instance
//MembraneAbinitioOP
moves::MoverOP
MembraneAbinitio::clone() const
{
	return moves::MoverOP( new MembraneAbinitio( *this ) );
}

void MembraneAbinitio::apply( pose::Pose & pose ) {
	using namespace moves;
	using namespace scoring;
	//bool success =
	Parent::apply( pose );
	total_trials_ = 0;

#ifdef GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( pose.conformation(), "start_pose" ); //add viewer
#endif

	//std::exit(1);
	if ( !only_stage4_ ) {
		if ( !recover_checkpoint( pose, "stage_1") ) {

			PROF_START( basic::STAGE1 );
			clock_t starttime = clock();
			// part 1 ----------------------------------------
			tr.Info <<  "\n===================================================================\n";
			tr.Info <<  "   Stage 1                                                         \n";
			tr.Info <<  "   Folding with score0 for max of " << stage1_cycles() << std::endl;

			prepare_stage1( pose );

			if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
				prof_show();
				output_debug_structure( pose, "stage0" );
			}
			//std::cout << "MOVABLE JUMP " << movemap_->get_jump(1);
			//std::cout << "\n";
			/*
			for(Size i=1;i<=pose.total_residue();++i) {
			pose.set_phi(i,-60);
			pose.set_psi(i,-40);
			if(movemap_->get_bb(i)){
			std::cout << 1;
			} else {
			std::cout << 0;
			}
			}
			std::cout <<"\n";
			*/
			//pose.dump_pdb("score0pre.pdb");
			//std::exit(1);
			do_stage1_cycles( pose );
			//pose.dump_pdb("score0.pdb");
			//std::exit(1);
			// mc().recover_low( pose ); seems to be a bad choice after score0

			if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();

			clock_t endtime = clock();
			PROF_STOP( basic::STAGE1 );

			if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
				tr << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
				prof_show();
				output_debug_structure( pose, "stage1" );
			}

			checkpoint( pose, "stage_1");
		}

		//  if ( option[ OptionKeys::abinitio::only_stage1_2 ] ) return success;


		// part 2 ----------------------------------------


		core::scoring::MembraneTopology const & topology(MembraneTopology_from_pose( pose ));
		nonconst_MembraneTopology_from_pose(pose).reset_allowed_scoring();
		Size total_tmhelix(topology.tmhelix());
		Size tmh_inserted(0);


		//ADD
		add_spanning_region(pose);
		while ( tmh_inserted<total_tmhelix )
				{
			tr.Info <<  "\n===================================================================\n";
			tr.Info <<  "   Stage 2                                                         \n";
			tr.Info <<  "   Starting sequencial insertion of " << total_tmhelix << ", membrane helices " << tmh_inserted << " inserted so far\n";
			tr.Info <<  "   Folding with score_membrane for " << stage2_cycles() << std::endl;
			print_debug(pose);
			{
				PROF_START( basic::STAGE2 );
				clock_t starttime = clock();


				prepare_stage2( pose );

				do_stage2_cycles( pose );
				mc().recover_low( pose );

				if  ( tr.visible() ) current_scorefxn().show( tr, pose );
				mc().show_counters();
				total_trials_+=mc().total_trials();
				mc().reset_counters();

				clock_t endtime = clock();
				PROF_STOP( basic::STAGE2 );

				if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
					// pose.dump_pdb("stage2.pdb");
					output_debug_structure( pose, "stage2" );
					tr << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
					prof_show();
				}

				//  checkpoint( pose, "stage_2" );
			}

			{
				//BW ADD ALLOW FRAGMENT INSERTION IN ALL INSERTED REGIONS
				move_all_inserted(pose);


				// moved checkpointing into do_stage3_cycles because of structure store
				// part 3 ----------------------------------------
				tr.Info <<  "\n===================================================================\n";
				tr.Info <<  "   Stage 3                                                         \n";
				tr.Info <<  "   Folding all inserted regions  with score_membrane for " << stage3_cycles() <<std::endl;
				//     tr.Info <<  std::endl << movemap_;
				print_debug(pose);
				PROF_START( basic::STAGE3 );
				clock_t starttime = clock();

				prepare_stage3( pose );
				if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );

				do_stage3_cycles( pose );
				mc().recover_low( pose );

				if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );
				mc().show_counters();
				total_trials_+=mc().total_trials();
				mc().reset_counters();

				clock_t endtime = clock();
				PROF_STOP( basic::STAGE3);

				if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
					output_debug_structure( pose, "stage3" );
					tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC) << std::endl;
					prof_show();
				}

				//add the extra score25 (stage3b)
				{
					// moved checkpointing into do_stage3_cycles because of structure store
					// part 3 ----------------------------------------
					tr.Info <<  "\n===================================================================\n";
					tr.Info <<  "   Stage 3b                                                         \n";
					tr.Info <<  "   Folding all inserted regions with score_membrane for " << stage3_cycles() <<std::endl;
					//     tr.Info <<  std::endl << movemap_;
					print_debug(pose);
					PROF_START( basic::STAGE3 );
					clock_t starttime = clock();

					prepare_stage3( pose );
					if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );

					do_stage3b_cycles( pose );
					mc().recover_low( pose );

					if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );
					mc().show_counters();
					total_trials_+=mc().total_trials();
					mc().reset_counters();

					clock_t endtime = clock();
					PROF_STOP( basic::STAGE3);

					if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
						output_debug_structure( pose, "stage3" );
						tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC) << std::endl;
						prof_show();
					}
				}

			}
			//ADD A NEW REGION TO BE INSERTED
			//tr.Info << "Before: tmh_inserted total_tmhelix " << tmh_inserted << " " << topology.tmh_inserted() << " " << total_tmhelix << std::endl;
			tmh_inserted=get_tmh_inserted(pose);
			//tr.Info << "Before: tmh_inserted total_tmhelix " << tmh_inserted << " " << topology.tmh_inserted() << " " << total_tmhelix << std::endl;
			add_spanning_region(pose);
			//tr.Info << "After: tmh_inserted total_tmhelix " << tmh_inserted << " " << topology.tmh_inserted() << " " << total_tmhelix << std::endl;
		}
	} // if (! only_stage4 )
	move_all_inserted(pose);


	// part 4 ------------------------------------------
	tr.Info <<  "\n===================================================================\n";
	tr.Info <<  "   Stage 4                                                         \n";
	tr.Info <<  "   Folding ALL REGIONS with score_membrane for " << stage4_cycles() <<std::endl;
	print_debug(pose);
	PROF_START( basic::STAGE4 );
	clock_t starttime = clock();

	prepare_stage4( pose );

	if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);

	do_stage4_cycles( pose );
	mc().recover_low( pose );

	if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);
	mc().show_counters();
	total_trials_+=mc().total_trials();
	mc().reset_counters();

	clock_t endtime = clock();
	PROF_STOP( basic::STAGE4 );

	if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
		output_debug_structure( pose, "stage4" );
		tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC ) << std::endl;
		prof_show();
	}

	tr.Info <<  "\n===================================================================\n";
	tr.Info <<  "   Finished Abinitio                                                 \n";
	tr.Info <<  "   Total trials " << total_trials_ << "\n";
	tr.Info <<  std::endl;

	// score_stage3a_->show( tr, pose );
	// score_stage3b_->show( tr, pose );
	// tr.flush();
	clear_checkpoints();
	// structure_store().clear();
	return; // true;
} // MembraneAbinitio::apply( pose::Pose & pose )


//@brief return FramgentMover for smooth_small fragment insertions (i.e., stage4 moves)
simple_moves::FragmentMoverOP
MembraneAbinitio::smooth_move_small() {
	return smooth_move_small_;
}

std::string
MembraneAbinitio::get_name() const {
	return "MembraneAbinitio";
}

//@brief return FragmentMover for small fragment insertions ( i.e., stage3/4 moves )
simple_moves::FragmentMoverOP
MembraneAbinitio::brute_move_small() {
	return brute_move_small_;
}

simple_moves::FragmentMoverOP
MembraneAbinitio::brute_move_small_top25() {
	return brute_move_small_top25_;
}

//@brief return FragmentMover for large fragment insertions (i.e., stage1/2 moves )
simple_moves::FragmentMoverOP
MembraneAbinitio::brute_move_large() {
	return brute_move_large_;
}

//@brief change the movemap ( is propagated to mover-objects )
//@detail overload if your extension stores additional moves as member variables
void
MembraneAbinitio::set_movemap( core::kinematics::MoveMapCOP mm )
{
	movemap_ = mm;
	if ( smooth_move_small_ ) smooth_move_small_ -> set_movemap( mm );
	if ( brute_move_small_  ) brute_move_small_  -> set_movemap( mm );
	if ( brute_move_small_top25_  ) brute_move_small_top25_  -> set_movemap( mm );
	if ( brute_move_large_  ) brute_move_large_  -> set_movemap( mm );
}

//@brief set new instances of FragmentMovers
void
MembraneAbinitio::set_moves(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small
)
{
	smooth_move_small_ = smooth_move_small;
	brute_move_small_  = brute_move_small;
	brute_move_large_  = brute_move_large;
	update_moves();
}

//@brief returns current movemap
core::kinematics::MoveMapCOP
MembraneAbinitio::movemap() {
	return movemap_;
}

//@detail read cmd_line options and set default versions for many protocol members: trials/moves, score-functions, Monte-Carlo
void MembraneAbinitio::set_defaults( pose::Pose const& pose ) {
	temperature_ = 2.0;
	set_default_scores();
	set_default_options();
	set_default_mc( pose, *score_stage1_ );
	update_moves();
}

//@detail called to notify about changes in Movers: new movemap or Moverclass
void MembraneAbinitio::update_moves() {
	/* set apply_large_frags_ and
	short_insert_region_
	*/
	/* what about move-map ? It can be set manually for all Fragment_Moves .. */
	// set_move_map();
	set_trials();
}

//@detail create instances of TrialMover for our FragmentMover objects
void MembraneAbinitio::set_trials() {
	// setup loop1
	assert( brute_move_large_ );
	trial_large_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_large_, mc_ ) );
	// trial_large_->keep_stats_type( false );
	//trial_large_->keep_stats_type( moves::all_stats );

	assert( brute_move_small_ );
	trial_small_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_small_, mc_ ) );
	//trial_small_->set_keep_stats( false );
	// trial_small_->keep_stats_type( false );
	//trial_small_->keep_stats_type( moves::all_stats );

	assert( brute_move_small_top25_ );
	trial_small_top25_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_small_top25_, mc_ ) );
	//trial_small_top25_ ->set_keep_stats( false );
	//trial_small_top25_ ->keep_stats_type( false );
	//trial_small_top25_ ->keep_stats_type( moves::all_stats );

	assert( smooth_move_small_ );
	smooth_trial_small_ = moves::TrialMoverOP( new moves::TrialMover( smooth_move_small_, mc_ ) );
	// smooth_trial_small_->set_keep_stats( false );
	//smooth_trial_small_->keep_stats_type( false );
	//smooth_trial_small_ ->keep_stats_type( moves::all_stats );
}

//@detail sets Monto-Carlo object to default
void MembraneAbinitio::set_default_mc( pose::Pose const& pose, scoring::ScoreFunction const& scorefxn ) {
	set_mc( moves::MonteCarloOP( new moves::MonteCarlo( pose, scorefxn, temperature_ ) ) );
}

//@detail sets Monto-Carlo object
void MembraneAbinitio::set_mc( moves::MonteCarloOP mc_in ) {
	mc_ = mc_in;
	if ( trial_large_ ) trial_large_->set_mc( mc_ );
	if ( trial_small_ ) trial_small_->set_mc( mc_ );
	if ( trial_small_top25_ ) trial_small_top25_ ->set_mc( mc_ );
	if ( smooth_trial_small_ ) smooth_trial_small_->set_mc( mc_ );
}

//@detail override cmd-line setting for "increase_cycling"
void MembraneAbinitio::set_cycles( Real increase_cycles ) {
	stage1_cycles_ = static_cast< int > (10000 * increase_cycles); //it will bail out if all fragments are replaced and since we are dealing with proteins longer than 50 residues we need some more cycles.
	stage2_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage3_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage4_cycles_ = static_cast< int > (4000 * increase_cycles);

	using namespace basic::options;
	// if ( option[ OptionKeys::abinitio::only_stage1_2 ]() ) {
	//  stage2_cycles_ = 0;
	//  stage3_cycles_ = 0;
	//  stage4_cycles_ = 0;
	// }
}

void MembraneAbinitio::set_default_scores() {
	using namespace scoring;
	tr.Debug << "creating membrane scoring functions" << std::endl;
	score_stage1_  = ScoreFunctionFactory::create_score_function( "score0_membrane" );
	score_stage2_  = ScoreFunctionFactory::create_score_function( "score_membrane" );
	//score_stage2_  = ScoreFunctionFactory::create_score_function( "score2" );
	score_stage3a_ = ScoreFunctionFactory::create_score_function( "score_membrane" );
	score_stage3b_ = ScoreFunctionFactory::create_score_function( "score_membrane" );
	score_stage4_  = ScoreFunctionFactory::create_score_function( "score_membrane" );

}

/// @brief sets a score weight for all stages of abinitio
void MembraneAbinitio::set_score_weight( scoring::ScoreType type, Real setting, StageID stage ) {
	tr.Debug << "set score weights for ";
	if ( stage == ALL_STAGES ) tr.Debug << "all stages ";
	else tr.Debug << "stage " << (stage <= STAGE_3a ? stage : ( stage-1 ) ) << ( stage == STAGE_3b ? "b " : " " );
	tr.Debug << scoring::name_from_score_type(type) << " " << setting << std::endl;
	if ( score_stage1_  && ( stage == STAGE_1  || stage == ALL_STAGES ) ) score_stage1_ ->set_weight(type, setting);
	if ( score_stage2_  && ( stage == STAGE_2  || stage == ALL_STAGES ) ) score_stage2_ ->set_weight(type, setting);
	if ( score_stage3a_ && ( stage == STAGE_3a || stage == ALL_STAGES ) ) score_stage3a_->set_weight(type, setting);
	if ( score_stage3b_ && ( stage == STAGE_3b || stage == ALL_STAGES ) ) score_stage3b_->set_weight(type, setting);
	if ( score_stage4_  && ( stage == STAGE_4  || stage == ALL_STAGES ) ) score_stage4_ ->set_weight(type, setting);
}

//@brief currently used score function ( depends on stage )
scoring::ScoreFunction const& MembraneAbinitio::current_scorefxn() const {
	return mc().score_function();
}

//@brief set current scorefunction
void MembraneAbinitio::current_scorefxn( scoring::ScoreFunction const& scorefxn ) {
	mc().score_function( scorefxn );
}

//@brief set individual weight of current scorefunction --- does not change the predifined scores: score_stageX_
void MembraneAbinitio::set_current_weight( core::scoring::ScoreType type, core::Real setting ) {
	scoring::ScoreFunctionOP scorefxn ( mc().score_function().clone() );
	scorefxn->set_weight( type, setting );
	mc().score_function( *scorefxn ); //trigger rescore
}

void MembraneAbinitio::set_default_options() {
	using namespace basic::options;
	just_smooth_cycles_ = option[ OptionKeys::abinitio::smooth_cycles_only ]; // defaults to false
	//  bQuickTest_ = option[ OptionKeys::abinitio::test_2 ];
	bQuickTest_ = false;
	if ( bQuickTest() ) {
		set_cycles( 0.001 );
	} else {
		set_cycles( option[ OptionKeys::abinitio::increase_cycles ] ); // defaults to factor of 1.0
	}

	only_stage4_ = just_smooth_cycles_ ;

	//bw fix for uninitialized value
	apply_large_frags_   = true;  // apply large frags in phase 2!

	// in rosetta++ switched on in fold_abinitio if contig_size < 30 in pose_abinitio never
	//bw fix for uninitialized value
	short_insert_region_ = false;  // apply small fragments in phase 2!
}


void MembraneAbinitio::add_spanning_region( core::pose::Pose & pose) {

	//using namespace common_regions;
	//using namespace misc;
	//using namespace membrane;
	//using namespace protein_maps;

	//bool exist;
	//bool none_selected = true;
	Size new_membrane_region=0;
	bool done=false;
	Size TMHs=0;
	//Size membrane_jump_counter=1;
	Size const num_cutpoints(pose.fold_tree().num_cutpoint());
	Size const num_jumps(pose.num_jump());
	core::scoring::MembraneTopology & topology( nonconst_MembraneTopology_from_pose(pose) );
	Size const total_tmhelix=topology.tmhelix();
	Size const total_residue=pose.total_residue();


	tr.Info << "Adding region tmh:" << total_tmhelix << " tmh_inserted: " << topology.tmh_inserted() << " nres: " << total_residue << " num_jumps: " << num_jumps << std::endl;

	//FArray2D_int const & jump_point( pose.fold_tree().get_jump_point() );

	utility::vector1< Size > const & cuts( pose.fold_tree().cutpoints() ); //( num_cutpoints ));

	//bw just define this vector to interface with the r++ code might change later...
	FArray2D_int jump_point(num_cutpoints,2);


	FArray1D_bool new_region(total_residue,false);
	FArray1D_bool inserted_regions(total_residue,false);

	//Assigns a vector with the helices numbered.
	FArray1D<Size> res_TMH_mapping(pose.total_residue()); //should perhaps be global..
	FArray2D_bool nb_tmh_list(total_tmhelix,total_tmhelix,false);

	//check all helices that are connect by a jump.
	FArray2D_bool tmh_jump(total_tmhelix,total_tmhelix,false);
	for ( Size j = 1; j <= pose.total_residue(); ++j ) {
		new_region(j) = false;
		inserted_regions(j) = false;

		//bw change definition of membrane region to include jumps to non-tmh.
		if ( j<=topology.span_end(1) ) { //membrane_helix(1,2))
			res_TMH_mapping(j)=1;
		} else if ( j>topology.span_end(total_tmhelix) ) {
			res_TMH_mapping(j)=total_tmhelix;
		} else {
			for ( Size reg = 2; reg <= total_tmhelix; ++reg ) {
				if ( j>topology.span_end(reg-1) && j<=topology.span_end(reg) ) { //membrane_helix( reg-1, 2 ) && j<=membrane_helix(reg,2))
					res_TMH_mapping(j)=reg;
				}
			}
		}
		//  tr.Info << "RES TMH " << j << " " << res_TMH_mapping(j) << std::endl;
	}

	//mark tmh as neighbours if they are consecutive and no cut is between them
	for ( Size reg = 1; reg < total_tmhelix; ++reg ) {
		// check cuts.
		bool no_cut=true;
		for ( Size i=topology.span_end(reg); i<topology.span_begin(reg+1); ++i ) { // (int i=membrane_helix(reg,2);i<membrane_helix(reg+1,2)??? no cut in tmh so it should be the same as check to begin;++i){
			for ( Size j=1; j<=num_cutpoints; ++j ) {
				if ( cuts[j]==i ) {
					no_cut=false;
				}
			}
		}
		if ( no_cut ) {
			nb_tmh_list(reg,reg+1)=true;
			nb_tmh_list(reg+1,reg)=true;
			//std::cout << "NO  " << reg << "," << reg+1 << std::endl;
		}
	}


	//check all helices that are connect by a jump.
	//diagonal elements are true if the helix is involved in a jump.
	for ( Size i=1; i<=num_jumps; ++i ) {
		Size d=pose.fold_tree().downstream_jump_residue(i);
		Size u=pose.fold_tree().upstream_jump_residue(i);
		tmh_jump(res_TMH_mapping(d),res_TMH_mapping(u))=true;
		tmh_jump(res_TMH_mapping(d),res_TMH_mapping(d))=true;
		tmh_jump(res_TMH_mapping(u),res_TMH_mapping(d))=true;
		tmh_jump(res_TMH_mapping(u),res_TMH_mapping(u))=true;
		nb_tmh_list(res_TMH_mapping(d),res_TMH_mapping(u))=true;
		nb_tmh_list(res_TMH_mapping(u),res_TMH_mapping(d))=true;

	}
	// check if all helices involved in a jump have been inserted.
	// produce a vector with insertable regions.
	// a region is insertable if it is next to something that is already inserted.
	// either by sequence or by a jump. Or if nothing is inserted everything is insertable.

	// Choose to insert a jump region first. The next jump to be inserted will be choosen randomly.

	// When the fold tree gets complicted adjecent helices might not be in contact.

	FArray1D<Size> insertable_region(total_tmhelix);
	Size k=0;
	if ( topology.tmh_inserted()==0 ) {
		for ( Size i = 1; i <= total_tmhelix; ++i ) {
			if ( tmh_jump(i,i) && !topology.allow_tmh_scoring(i) ) { //!TMH_done(i))
				++k;
				insertable_region(k)=i;
				tr.Info << "INSERTABLE JUMP REGION: " << i << std::endl;
			}
		}
	}
	if ( k==0 ) { // if no jumps are available
		//non jump helices
		if ( topology.tmh_inserted()==0 ) { //first pass?
			for ( Size j = 1; j <= total_tmhelix; ++j ) {
				++k;
				insertable_region(k)=j;
				tr.Info << "INSERTABLE REGION FIRST PASS: " << j << std::endl;
			}
		} else {
			for ( Size i = 1; i <= total_tmhelix; ++i ) {
				if ( topology.allow_tmh_scoring(i) ) {
					for ( Size j = 1; j <= total_tmhelix; ++j ) {
						if ( !topology.allow_tmh_scoring(j) && nb_tmh_list(i,j) ) { // neigbor in fold tree
							++k;
							insertable_region(k)=j;
							tr.Info << "INSERTABLE REGION: " << j << std::endl;
						}
					}
				}
			}
		}
	}
	if ( k>0 ) { // any more regions to insert?
		Size index=1 + static_cast< int >( numeric::random::uniform()*k);
		FArray1D_bool new_membrane_regions(total_tmhelix,false);  //need an array to handle the case when more than helix is inserted at the time.
		new_membrane_region=insertable_region(index);
		new_membrane_regions(new_membrane_region)=true;
		if ( tmh_jump(new_membrane_region,new_membrane_region) ) { // a "jump helix"
			//insert all helices that are connected through jumps
			// could in principle be a whole network of jumps..
			for ( Size i=1; i<=total_tmhelix; ++i ) {
				new_membrane_regions(i)=tmh_jump(new_membrane_region,i);
				if ( new_membrane_regions(i) ) {
					for ( Size j=1; j<=total_tmhelix; ++j ) {
						new_membrane_regions(j)=tmh_jump(i,j);
					}
				}
			}
			for ( Size i=1; i<=total_tmhelix; ++i ) {
				if ( new_membrane_regions(i) ) {
					topology.set_allow_tmh_scoring(i,true);
					tr.Info << "NEW JUMP REGION: index " << i << std::endl;
				}
			}
		} else {

			if ( topology.tmh_inserted()==0 ) { // FIRST PASS
				Size new_membrane_region2;
				if ( new_membrane_region==1 ) {
					new_membrane_region2=new_membrane_region+1;
				} else if ( new_membrane_region==total_tmhelix ) {
					new_membrane_region2=new_membrane_region-1;
				} else if ( numeric::random::uniform() < 0.5 ) {
					new_membrane_region2=new_membrane_region+1;
				} else {
					new_membrane_region2=new_membrane_region-1;
				}
				tr.Info << "NEW REGION FIRST PASS: region1 region2 " << new_membrane_region << " " << new_membrane_region2 << std::endl;
				new_membrane_regions(new_membrane_region2)=true;
				topology.set_allow_tmh_scoring(new_membrane_region2,true);

			}


			tr.Info << "NEW REGION: index " << new_membrane_region << " " << index << " " << k << std::endl;
			//   TMH_done(new_membrane_region)=true;
			topology.set_allow_tmh_scoring(new_membrane_region,true);
			//   allow_tmh_scoring( new_membrane_region ) = true;
			// define maximum range middle_helix_num1 should be the first and middle_helix_num2 should be the last.
			// membrane_score_quick then uses the the allow_tmh_scoring vector to know which helix to include in the
			// membrane normal calculation. Assume they are correctly set before hand. We know that they span the min-max range
			// and we only have to check if the newly added region will change the min-max boundery.
			//  if(new_membrane_region>middle_helix_num2)
			//   middle_helix_num2=new_membrane_region;
			//  if(new_membrane_region<middle_helix_num1)
			//   middle_helix_num1=new_membrane_region;
		}

		//Mark up which residues that are allowed to scored based on the TMH_done.
		Size start,end;
		TMHs=0;
		for ( Size i = 1; i <= total_tmhelix; ++i ) {
			if ( topology.allow_tmh_scoring(i) ) { //TMH_done(i)){
				if ( i==1 ) {
					start=1;
				} else {
					start=topology.span_end(i-1); //membrane_helix( i-1, 2 ); //from old code
				}
				if ( i==total_tmhelix ) {
					end=pose.total_residue();
				} else {
					end=topology.span_begin(i+1); //membrane_helix( i+1, 1 ); //from old code
				}
				for ( Size j=1; j<=num_cutpoints; ++j ) {
					if ( cuts[j] >= start &&  cuts[j] <= end ) {
						if ( cuts[j]-start >= end-cuts[j] ) {
							end=cuts[j]; //the residue is cut off so do not include it in scoring.
						} else {
							start=cuts[j]+1;
						}
					}
				}
				for ( Size j=start; j<=end; ++j ) {
					inserted_regions(j)=true;
					if ( new_membrane_regions(i) ) {
						new_region(j)=true;
					}
				}
				++TMHs;
			}
		}
		done=false;
	} else {
		done=true;
		for ( Size j = 1; j <= total_residue; ++j ) {
			inserted_regions(j) = true;
		}
		TMHs=total_tmhelix;
	}
	topology.set_tmh_inserted(TMHs);
	kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
	if ( done ) {
		new_mm->set_bb( true );
		for ( Size i = 1; i <= total_residue; ++i ) {
			topology.set_allow_scoring(i,true);
		}
	} else {
		for ( Size i = 1; i <= total_residue; ++i ) {
			if ( inserted_regions(i) ) {
				topology.set_allow_scoring(i,true);
			} else {
				topology.set_allow_scoring(i,false);
			}
			if ( new_region(i) ) {
				new_mm->set_bb(i,true);
			} else {
				new_mm->set_bb(i,false);
			}
		}
	}
	set_movemap(new_mm);

	/*
	middle_helix_num1=1;
	middle_helix_num2=0;
	for ( int j = 1; j <= total_tmhelix; ++j ) {
	allow_tmh_scoring(j)=false;
	if(TMH_done(j)) {
	allow_tmh_scoring(j)=true;
	middle_helix_num2=j;
	} else if(middle_helix_num2==0) {
	middle_helix_num1=j;
	}
	}
	//bw This are important for scoring still...


	mres_start=1;
	mres_end=1;
	for ( int j = 1; j <= pose.total_residue(); ++j ) {
	if(inserted_regions(j))
	{
	mres_end=j;
	}
	else if(mres_end==1)
	{
	mres_start=j;
	}
	}

	//  exist = reset_insert_map();
	//for ( int j = 1; j <= total_residue; ++j ) {
	//  pose.set_allow_bb_move(j,allow_insert(j));
	//}

	*/
	tr.Info << "max_tm new_reg TMHs done tmh_inserted " << std::endl;

	tr.Info << I( 7, total_tmhelix ) << ' ' << I( 7, new_membrane_region ) << ' ' << I(4, TMHs) << ' ' << done << ' ' << topology.tmh_inserted() << std::endl;
	/*
	<< ' ' << I( 7, middle_helix_num1 )
	<< ' ' << I( 7, middle_helix_num2 ) << ' ' << I( 7, mres_start ) << ' ' << I( 8, mres_end ) << ' '
	<< I( 4, TMHs ) << ' ' << L( 4, done ) << "  " << A( 6, TMHpred_method ) << std::endl;
	*/
	//for(int a=1;a<total_residue;++a)
	// {
	//  if(allow_insert(a))
	//   {
	//    std::cout << "1";
	//   }
	//  else
	//   {
	//    std::cout << "0";
	//   }
	// }
	// std::cout << std::endl;
	//std::cout << insert_map(total_insert) << " " << total_insert << std::endl;


}

void MembraneAbinitio::move_all_inserted( core::pose::Pose & pose) {
	core::scoring::MembraneTopology & topology( nonconst_MembraneTopology_from_pose( pose ));
	kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( topology.allow_scoring(i) ) {
			new_mm->set_bb(i,true);
		} else {
			new_mm->set_bb(i,false);
		}
	}
	set_movemap(new_mm);
}

void MembraneAbinitio::print_debug(core::pose::Pose & pose)
{
	//  return;
	//  if(!option[basic::options::OptionKeys::abinitio::membrane_print])
	//   return;
	core::scoring::MembraneTopology const & topology( MembraneTopology_from_pose( pose ));
	Size total_tmhelix=topology.tmhelix();
	FArray1D_int res_TMH_mapping(pose.total_residue()); //should perhaps be global..
	for ( Size j = 1; j <= pose.total_residue(); ++j ) {
		Size tmh=0;
		for ( Size reg = 1; reg <= total_tmhelix; ++reg ) {
			if ( j>=topology.span_begin(reg) && j<=topology.span_end(reg) ) { //membrane_helix( reg-1, 2 ) && j<=membrane_helix(reg,2))
				tmh=reg;

			}
		}

		res_TMH_mapping(j)=tmh;
		//  tr.Info << "RES TMH " << j << " " << res_TMH_mapping(j) << std::endl;
	}

	//  pose.dump_pdb("debug.pdb");
	std::cout << "TMH     : ";


	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		std::cout << res_TMH_mapping(i);
	}
	std::cout <<"\n";
	std::cout << "SCORING : ";
	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( topology.allow_scoring(i) ) {
			std::cout << 1;
		} else {
			std::cout << 0;
		}
	}
	std::cout <<"\n";
	std::cout << "MOVEMAP : ";
	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( movemap_->get_bb(i) ) {
			std::cout << 1;
		} else {
			std::cout << 0;
		}
	}
	std::cout <<"\n";
	mc_->show_scores();
	current_scorefxn().show( std::cout, pose );

}

/// @brief (helper) functor class which keeps track of initial phi/psi values.
/// @detail
/// calls of operator ( pose ) compare the initial phi/psi values
////to the current values of the given pose. Returns false once all phi/psi values
/// have been modified.
class AllResiduesChanged : public moves::PoseCondition {
public:
	AllResiduesChanged( core::pose::Pose const & pose, core::fragment::InsertMap const& insert_map ) :
		insert_pos_( pose.total_residue(), false )
	{
		set_initial_pose( pose );
		compute_insert_pos( insert_map );
	}
private:

	void compute_insert_pos( core::fragment::InsertMap const& insert_map ) {
		for ( core::fragment::InsertMap::const_iterator it = insert_map.begin(),
				eit = insert_map.end(); it != eit; ++it ) {
			Size const pos ( *it );
			if ( pos > insert_pos_.size() ) break;
			insert_pos_[ pos ] = true;
		}
	}

	void set_initial_pose( const core::pose::Pose & pose ) {
		for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
			initial_phis.push_back( pose.phi(i) );
			initial_psis.push_back( pose.psi(i) );
		}

		original_sequence_ = pose.sequence();
	}

public:
	virtual bool operator() ( const core::pose::Pose & pose ) {
		assert( original_sequence_ == pose.sequence() ); // imperfect attempt to check that Pose hasn't changed ...
		for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
			if ( initial_phis[i] == pose.phi(i) && insert_pos_[ i ] ) {
				return false;
			}
			if ( initial_psis[i] == pose.psi(i) && insert_pos_[ i ] ) {
				return false;
			}
		}
		return true;
	}

private:
	utility::vector1< core::Real > initial_phis;
	utility::vector1< core::Real > initial_psis;

	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool initialized_;

	std::string original_sequence_;
	utility::vector1< bool > insert_pos_;
};

/*
/// @brief (helper) functor class which keeps track of old pose for the
/// convergence check in stage3 cycles
/// @detail
/// calls of operator ( pose ) compare the
class MonteCarloExceptionConverge;
typedef  utility::pointer::owning_ptr< MonteCarloExceptionConverge >  MonteCarloExceptionConvergeOP;

class MonteCarloExceptionConverge : public moves::PoseCondition {
public:
MonteCarloExceptionConverge() : bInit_( false ), ct_( 0 ) {};
void reset() { ct_ = 0; bInit_ = false; };
void set_trials( moves::TrialMoverOP trin ) {
trials_ = trin;
assert( trials_->keep_stats() );
last_move_ = 0;
};
virtual bool operator() ( const core::pose::Pose & pose );
private:
pose::Pose very_old_pose_;
bool bInit_;
Size ct_;
moves::TrialMoverOP trials_;
Size last_move_;
};

// keep going --> return true
bool MonteCarloExceptionConverge::operator() ( const core::pose::Pose & pose ) {
if ( !bInit_ ) {
bInit_ = true;
very_old_pose_ = pose;
return true;
}
assert( trials_ );
tr.Trace << "TrialCounter in MonteCarloExceptionConverge: " << trials_->num_accepts() << std::endl;
if ( numeric::mod(trials_->num_accepts(),100) != 0 ) return true;
if ( (Size) trials_->num_accepts() <= last_move_ ) return true;
last_move_ = trials_->num_accepts();
// change this later to this: (after we compared with rosetta++ and are happy)
// if ( numeric::mod(++ct_, 1000) != 0 ) return false; //assumes an approx acceptance rate of 0.1

// still here? do the check:

core::Real converge_rms = core::scoring::CA_rmsd( very_old_pose_, pose );
very_old_pose_ = pose;
if ( converge_rms >= 3.0 ) {
return true;
};
// if we get here thing is converged stop the While-Loop
tr.Info << " stop cycles in stage3 due to convergence " << std::endl;
return false;
}

*/
int MembraneAbinitio::do_stage1_cycles( pose::Pose &pose ) {
	AllResiduesChanged done( pose, brute_move_large()->insert_map() );
	Size j;
	for ( j = 1; j <= stage1_cycles(); ++j ) {
		// std::cout << j << " " << pose.energies() << "\n";

		trial_large()->apply( pose );// apply a large fragment insertion, accept with MC boltzmann probability
		//   bool tmp=done(pose);
		//   std::cout << tmp << "\n";
		//pose::Pose tmp_pose;

		//pose.dump_pdb(string_of(j)+".pdb");
		if ( j % 1000==0 ) {
			tr.Info << "Step (score0) " << j << std::endl;
		}
		if ( done(pose) ) {
			tr.Info << "Replaced extended chain after " << j << " cycles" << std::endl;
			mc().reset( pose ); // make sure that we keep the final structure
			return j;
		}
	}

	tr.Info << "Warning: extended chain may still remain after " << stage1_cycles() << " cycles!" << std::endl;
	mc().reset( pose ); // make sure that we keep the final structure
	return j;
}

int MembraneAbinitio::do_stage2_cycles( pose::Pose &pose ) {

	//setup cycle

	//if ( apply_large_frags_   ) cycle->add_mover( trial_large_->mover() );
	//if ( short_insert_region_ ) cycle->add_mover( trial_small_->mover() );

	/*
	Size j;
	Size cycles=stage2_cycles();
	//Size tmh_inserted=get_tmh_inserted(pose);
	// if(get_tmh_inserted(pose)==1) {
	//  cycles*=;
	// }

	std::cout << "Start stage2\n";
	for ( j = 1; j <= cycles; ++j ) {
	trial_large()->apply( pose );

	//if(option[basic::options::OptionKeys::abinitio::membrane_print] && j%100 == 0) {


	// std::cout << "large j " << j <<"\n";
	// print_debug(pose);
	//}

	//if( tmh_inserted > 1)
	// pose.dump_pdb("stage_2_TMHs"+ObjexxFCL::string_of( tmh_inserted )+"-"+ObjexxFCL::string_of( j )); //+"_"+ObjexxFCL::string_of(lct2)
	}
	std::cout << "Start 3mer\n";
	for ( j = 1; j <= cycles; ++j ) {
	trial_small_top25()->apply( pose );

	if(option[basic::options::OptionKeys::abinitio::membrane_print] && j%100 == 0) {
	std::cout << "small1 j " << j <<"\n";
	print_debug(pose);
	}


	}
	std::cout << "Start 3mer\n";
	for ( j = 1; j <= cycles; ++j ) {

	trial_small_top25()->apply( pose );

	if(option[basic::options::OptionKeys::abinitio::membrane_print] && j%100 == 0) {
	std::cout << "small2 j " << j <<"\n";
	print_debug(pose);
	}

	}
	clock_t stop = clock();
	double t = (double) (stop-start)/CLOCKS_PER_SEC;
	std::cout << "End stage2: time " << t << "\n";
	*/
	// clock_t start = clock();
	moves::SequenceMoverOP cycle( new moves::SequenceMover() );
	cycle->add_mover( trial_large_->mover() );
	cycle->add_mover( trial_small_top25_->mover() );
	cycle->add_mover( trial_small_top25_->mover() );

	Size nr_cycles = stage2_cycles(); // ( short_insert_region ? 2 : 1 );
	moves::TrialMoverOP trials( new moves::TrialMover( cycle, mc_ptr() ) );
	moves::RepeatMover( trials, nr_cycles ).apply(pose);
	// clock_t stop = clock();
	// double t = (double) (stop-start)/CLOCKS_PER_SEC;
	// std::cout << "End stage2: time " << t << "\n";
	// std::exit(1);
	// moves::RepeatMover( stage2_mover( pose, trials ), nr_cycles ).apply(pose);
	//*/
	//is there a better way to find out how many steps ? for instance how many calls to scoring?
	return 3*nr_cycles; //stage2_cycles(); //nr_cycles; // as best guess
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

*/int MembraneAbinitio::do_stage3_cycles( pose::Pose &pose ) {
	// interlaced score2 / score 5 loops

	// nloops1 and nloops2 could become member-variables and thus changeable from the outside
	int nloop1 = 1;
	int nloop2 = 9; //10; //careful: if you change these the number of structures in the structure store changes.. problem with checkpointing
	// individual checkpoints for each stage3 iteration would be a remedy. ...
	/*
	if ( short_insert_region_ ) {
	nloop1 = 2;
	nloop2 = 5;
	}

	MonteCarloExceptionConvergeOP convergence_checker ( NULL );
	if ( !option[ basic::options::OptionKeys::abinitio::skip_convergence_check ] ) {
	convergence_checker = new MonteCarloExceptionConverge;
	}
	*/
	Size cycles=stage3_cycles();
	//  if(get_tmh_inserted(pose)==1) {
	//   cycles*=0.01;
	//  }
	moves::TrialMoverOP trials = trial_large();
	int iteration = 1;
	for ( int lct1 = 1; lct1 <= nloop1; lct1++ ) {
		// if ( lct1 > 1 ) trials = trial_small(); //only with short_insert_region!
		for ( int lct2 = 1; lct2 <= nloop2; lct2++, iteration++  ) {
			tr.Debug << "Loop: " << lct1 << "   " << lct2 << std::endl;
			if ( lct2>3 ) trials = trial_small_top25();
			if ( !recover_checkpoint( pose, "stage_3_iter"+ObjexxFCL::string_of( lct1)+"_"+ObjexxFCL::string_of(lct2) ) ) {
				prepare_loop_in_stage3( pose, iteration, nloop1*nloop2 );// bw only sets temperature currently...

				// interlace score2/score5
				if ( numeric::mod( lct2, 2 ) == 0 || lct2 > 7 ) { // chainbreak ramping...
					Real chainbreak_weight=lct1*lct2*0.25;
					mc_->recover_low( pose );
					tr.Debug << "select score_stage3a..." << std::endl;
					( *score_stage3a_).set_weight(scoring::linear_chainbreak,chainbreak_weight);
					mc_->score_function( *score_stage3a_ );
					// mc_->reset( pose );
				} else {
					Real chainbreak_weight=lct1*lct2*0.05;
					mc_->recover_low( pose );
					tr.Debug << "select score_stage3b..." << std::endl;
					( *score_stage3b_).set_weight(scoring::linear_chainbreak,chainbreak_weight);
					mc_->score_function( *score_stage3b_ );
					// mc_->reset( pose );
				}

				tr.Debug << "  Score stage3 loop iteration " << lct1 << " " << lct2 << std::endl;
				/*  if ( convergence_checker ) {
				moves::TrialMoverOP stage3_trials = stage3_mover( pose, lct1, lct2, trials );
				convergence_checker->set_trials( stage3_trials ); //can be removed late
				moves::WhileMover( stage3_trials, stage3_cycles(), convergence_checker ).apply( pose );
				} else {    //no convergence check -> no WhileMover */
				moves::RepeatMover( stage3_mover( pose, lct1, lct2, trials ), stage3_cycles() ).apply( pose );
				//moves::RepeatMover( trials, stage3_cycles() ).apply( pose );


				//     for(Size j=1;j<=cycles;++j) {
				//        trials->apply(pose);
				/*
				if(option[basic::options::OptionKeys::abinitio::membrane_print] && j%100 == 0) {
				std::cout << "3 j " << j << " " << lct2 << "\n";
				print_debug(pose);
				}
				*/
				//      }

				// }
				checkpoint( pose, "stage_3" );
				//     structure_store().push_back( mc_->lowest_score_pose() );
			}//recover_checkpoint
		}; // loop 2
	}; // loop 1
	return nloop1*nloop2*cycles; //stage3_cycles();
}


int MembraneAbinitio::do_stage3b_cycles( pose::Pose &pose ) {
	// interlaced score2 / score 5 loops

	// nloops1 and nloops2 could become member-variables and thus changeable from the outside
	int nloop1 = 1;
	int nloop2 = 2; //10; //careful: if you change these the number of structures in the structure store changes.. problem with checkpointing
	// individual checkpoints for each stage3 iteration would be a remedy. ...
	/*
	if ( short_insert_region_ ) {
	nloop1 = 2;
	nloop2 = 5;
	}

	MonteCarloExceptionConvergeOP convergence_checker ( NULL );
	if ( !option[ basic::options::OptionKeys::abinitio::skip_convergence_check ] ) {
	convergence_checker = new MonteCarloExceptionConverge;
	}
	*/Size cycles=stage3_cycles();
	//if(get_tmh_inserted(pose)==1) {
	// cycles*=0.01;
	//}
	moves::TrialMoverOP trials = trial_small_top25();
	int iteration = 1;
	for ( int lct1 = 1; lct1 <= nloop1; lct1++ ) {
		// if ( lct1 > 1 ) trials = trial_small(); //only with short_insert_region!
		for ( int lct2 = 1; lct2 <= nloop2; lct2++, iteration++  ) {
			tr.Debug << "Loop: " << lct1 << "   " << lct2 << std::endl;
			//   if(lct2>3) trial_small();
			if ( !recover_checkpoint( pose, "stage_3_iter"+ObjexxFCL::string_of( lct1)+"_"+ObjexxFCL::string_of(lct2) ) ) {
				prepare_loop_in_stage3( pose, iteration, nloop1*nloop2 );// bw only sets temperature currently...

				// interlace score2/score5
				if ( numeric::mod( lct2, 2 ) == 0 || lct2 > 7 ) { // chainbreak ramping...
					Real chainbreak_weight=lct1*lct2*0.25;
					mc_->recover_low( pose );
					tr.Debug << "select score_stage3a..." << std::endl;
					( *score_stage3a_).set_weight(scoring::linear_chainbreak,chainbreak_weight);
					mc_->score_function( *score_stage3a_ );
					// mc_->reset( pose );
				} else {
					Real chainbreak_weight=lct1*lct2*0.05;
					mc_->recover_low( pose );
					tr.Debug << "select score_stage3b..." << std::endl;
					( *score_stage3b_).set_weight(scoring::linear_chainbreak,chainbreak_weight);
					mc_->score_function( *score_stage3b_ );
					// mc_->reset( pose );
				}

				tr.Debug << "  Score stage3 loop iteration " << lct1 << " " << lct2 << std::endl;
				/*  if ( convergence_checker ) {
				moves::TrialMoverOP stage3_trials = stage3_mover( pose, lct1, lct2, trials );
				convergence_checker->set_trials( stage3_trials ); //can be removed late
				moves::WhileMover( stage3_trials, stage3_cycles(), convergence_checker ).apply( pose );
				} else {    //no convergence check -> no WhileMover */
				moves::RepeatMover( stage3_mover( pose, lct1, lct2, trials ), stage3_cycles() ).apply( pose );
				//moves::RepeatMover( trials, stage3_cycles() ).apply( pose );

				//for(Size j=1;j<=cycles;++j) {
				// trials->apply(pose);
				/*
				if(option[basic::options::OptionKeys::abinitio::membrane_print] && j%100 == 0) {

				std::cout << "3b j " << j << " " << lct2 << "\n";
				print_debug(pose);
				}*/

				//}
				// }
				checkpoint( pose, "stage_3" );
				//    structure_store().push_back( mc_->lowest_score_pose() );
			}//recover_checkpoint
		}; // loop 2
	}; // loop 1
	return nloop1*nloop2*cycles; //stage3_cycles();
}


moves::TrialMoverOP
MembraneAbinitio::stage2_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}


moves::TrialMoverOP
MembraneAbinitio::stage3_mover( pose::Pose &, int, int, moves::TrialMoverOP trials ) {
	return trials;
}

moves::TrialMoverOP
MembraneAbinitio::stage4_mover( pose::Pose &, int, moves::TrialMoverOP trials ) {
	return trials;
}

// interlaced score2 / score 5 loops
/*! @detail stage4 cycles:
nloop_stage4: iterations
stage4_cycle : trials per  iteration

first iteration: use trial_small()
following iterations: use trial_smooth()
only trial_smooth() if just_smooth_cycles==true

prepare_loop_in_stage4() is called each time before the stage4_cycles_ of trials are started.

start each iteration with the lowest_score_pose. ( mc->recover_low()  in prepare_loop_in_stage4()  )

*/
int MembraneAbinitio::do_stage4_cycles( pose::Pose &pose ) {
	Size nloop_stage4 = 3;

	for ( Size kk = 1; kk <= nloop_stage4; ++kk ) {
		if ( !recover_checkpoint( pose, "stage4_kk_" + ObjexxFCL::string_of(kk)) ) {
			//  if ( kk == 2 ) choose_frag_set_top_N_frags(number3merfrags);
			moves::TrialMoverOP trials;
			if ( kk == 1 && !just_smooth_cycles_ ) {
				trials = trial_small_top25();
			} else if ( kk==2 ) {
				trials = trial_small();
			} else {
				tr.Debug << "switch to smooth moves" << std::endl;
				trials = trial_smooth();
			}
			Real chainbreak_weight=kk*0.5+2.5;
			( *score_stage4_).set_weight(scoring::linear_chainbreak,chainbreak_weight);
			mc_->score_function( *score_stage4_ );
			mc_->reset( pose );
			tr.Debug << "prepare ..." << std::endl ;
			prepare_loop_in_stage4( pose, kk, nloop_stage4 );
			tr.Debug << "start " << stage4_cycles() << " cycles" << std::endl;
			moves::RepeatMover( stage4_mover( pose, kk, trials ), stage4_cycles() ).apply(pose);
			//don't store last structure since it will be exactly the same as the final structure delivered back via apply
			//   structure_store().push_back( mc_->lowest_score_pose() );
			tr.Debug << "finished" << std::endl;
			checkpoint( pose, "stage4_kk_" + ObjexxFCL::string_of(kk));
		}
	}  // loop kk
	return stage4_cycles() * nloop_stage4;
}

// anything you want to have done before the stages ?

// prepare stage1 sampling
void MembraneAbinitio::prepare_stage1( core::pose::Pose &pose ) {
	// don't do the following, because it destroys a user-defined mc.
	// set_default_mc( pose, *score_stage1_ ); //get a new mc-Object in case pose is not the same as before
	// Real chainbreak_weight=1.0;
	// ( *score_stage1_).set_weight(scoring::linear_chainbreak,chainbreak_weight);
	// ( *score_stage1_).set_weight(scoring::chainbreak,chainbreak_weight);

	mc_->score_function( *score_stage1_ );
	mc_->set_autotemp( false, temperature_ );
	mc_->set_temperature( temperature_ );
	mc_->reset( pose );
	(*score_stage1_)( pose );
	//score_stage1_->accumulate_residue_total_energies( pose ); // fix this
	kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
	new_mm->set_bb( true );
	set_movemap(new_mm);
}

void MembraneAbinitio::prepare_stage2( core::pose::Pose &pose ) {
	mc_->score_function( *score_stage2_ );
	mc_->set_autotemp( true, temperature_ );
	mc_->set_temperature( temperature_ ); // temperature might have changed due to autotemp..
	mc_->reset( pose );

	(*score_stage2_)(pose);
	///score_stage2_->accumulate_residue_total_energies( pose );
}


void MembraneAbinitio::prepare_stage3( core::pose::Pose &pose ) {
	mc_->score_function( *score_stage3a_ );
	mc_->set_autotemp( true, temperature_ );
	mc_->set_temperature( temperature_ ); // temperature might have changed due to autotemp..
	mc_->reset( pose );
	//score for this stage is changed in the do_stage3_cycles explicitly
	/*if ( option[ templates::change_movemap_2 ].user() && option[ templates::change_movemap_2 ] == 3 ) {
	kinematics::MoveMapOP new_mm = new kinematics::MoveMap( *movemap() );
	new_mm->set_bb( true );
	set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	*/
}


void MembraneAbinitio::prepare_stage4( core::pose::Pose &pose ) {
	mc_->set_autotemp( true, temperature_ );
	mc_->set_temperature( temperature_ ); // temperature might have changed due to autotemp..
	mc_->score_function( *score_stage4_ );
	mc_->reset( pose );

	(*score_stage4_)( pose );
	///score_stage4_->accumulate_residue_total_energies( pose ); // fix this
	/*
	if ( option[ templates::change_movemap_2 ].user() && option[ templates::change_movemap_2 ] == 4 ) {
	kinematics::MoveMapOP new_mm = new kinematics::MoveMap( *movemap() );
	new_mm->set_bb( true );
	set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	*/
}

void MembraneAbinitio::prepare_loop_in_stage3( core::pose::Pose &/*pose*/, Size, Size ){
	//  mc_->recover_low( pose );
	mc_->set_temperature( temperature_ );
}

void MembraneAbinitio::prepare_loop_in_stage4( core::pose::Pose &/*pose*/, Size, Size ){
	// mc_->recover_low( pose );
	mc_->set_temperature( temperature_ );
}

bool MembraneAbinitio::recover_checkpoint( pose::Pose &pose, std::string const & id ) {
	return checkpoint_.recover_checkpoint( pose, get_current_tag(), id );
}
//  if (!protocols::checkpoint::Timer::is_on()) return false; // required for checkpoint

//  // must be same as in checkpoint function
//  std::string checkpoint_id( type() + "_" + get_current_tag() + "_" + id );

//  if (utility::file::file_exists( checkpoint_id + ".pdb" )) {
//   utility::io::izstream izs(checkpoint_id + ".rng.state.gz");
//   numeric::random::RandomGenerator::restoreAllStates(izs);
//   izs.close();
//   core::import_pose::centroid_pose_from_pdb( pose, checkpoint_id + ".pdb" );
//   if (utility::file::file_exists( checkpoint_id + ".mc_last.pdb" )) {
//    pose::Pose recovered_mc_last;
//    core::import_pose::centroid_pose_from_pdb( recovered_mc_last, checkpoint_id + ".mc_last.pdb" );
//    mc_->set_last_accepted_pose( recovered_mc_last );
//   }
//   if (utility::file::file_exists( checkpoint_id + ".mc_low.pdb" )) {
//    pose::Pose recovered_mc_low;
//    core::import_pose::centroid_pose_from_pdb( recovered_mc_low, checkpoint_id + ".mc_low.pdb" );
//    mc_->set_lowest_score_pose( recovered_mc_low );
//   }
//   checkpoint_ids_.push_back( checkpoint_id );
//   //    tr << "Recovered checkpoint: " << checkpoint_id << std::endl;
//   return true;
//  }

//  return false;
//}

void MembraneAbinitio::checkpoint( pose::Pose &pose, std::string const & id ) {
	checkpoint_.checkpoint( pose, get_current_tag(), id );
}
//  if (!protocols::checkpoint::Timer::time_to_checkpoint()) return; // required for checkpoint

//  if (get_current_tag() == "NoTag") return;

//  std::string checkpoint_id( type() + "_" + get_current_tag() + "_" + id );

//  checkpoint_ids_.push_back( checkpoint_id );

//  utility::io::ozstream ozs(checkpoint_id + ".rng.state.gz");
//  numeric::random::RandomGenerator::saveAllStates(ozs);
//  ozs.close();

//  io::pdb::dump_pdb( mc_->last_accepted_pose(), checkpoint_id + ".mc_last.pdb" );
//  io::pdb::dump_pdb( mc_->lowest_score_pose(), checkpoint_id + ".mc_low.pdb" );
//  io::pdb::dump_pdb( pose, checkpoint_id + ".pdb" );
//  //  tr << "Created checkpoint: " << checkpoint_id << std::endl;

//  protocols::checkpoint::Timer::reset(); // required for checkpoint
// }

/// @brief for debugging, one wants to have access to the native pose.
//void MembraneAbinitio::set_native_pose( core::pose::Pose const & /*pose*/ ) {
// native_pose_ = pose; // this is bad since one doesn't know if set or not. make it a Pointer!
//}

void MembraneAbinitio::clear_checkpoints() {
	checkpoint_.clear_checkpoints();
}
//  for ( int i = 0; i < int( checkpoint_ids_.size() ); i++ ) {
//   //    tr << "deleting checkpoint files with id: " << checkpoint_ids_[i] << std::endl;
//   utility::file::file_delete( checkpoint_ids_[i] + ".mc_last.pdb" );
//   utility::file::file_delete( checkpoint_ids_[i] + ".mc_low.pdb" );
//   utility::file::file_delete( checkpoint_ids_[i] + ".pdb" );
//   utility::file::file_delete( checkpoint_ids_[i] + ".rng.state.gz" );
//  }
//  checkpoint_ids_.clear();
// } // MembraneAbinitio::clear_checkpoints()

void MembraneAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix ) {
	using namespace core::io::silent;
	if ( option[ basic::options::OptionKeys::out::file::silent ].user() ) {
		SilentFileData sfd;
		std::string silent_file = option[ basic::options::OptionKeys::out::file::silent ]() + "_" + prefix;

		mc().score_function()( pose );

		ProteinSilentStruct pss;
		pss.fill_struct( pose, get_current_tag() );

		evaluate_pose( pose, get_current_tag(), pss );

		sfd.write_silent_struct( pss, silent_file, !option[ OptionKeys::abinitio::debug_structures_2 ]()  /* bWriteScoresOnly */ );
	} // if option[ out::file::silent ].user()

	if ( option[ basic::options::OptionKeys::abinitio::explicit_pdb_debug ]() ) {
		pose.dump_pdb( prefix + get_current_tag() + ".pdb" );
	}
	/*
	if ( option[ basic::options::OptionKeys::abinitio::log_frags_2 ].user() ) {
	std::string filename = prefix + "_" + get_current_tag() + "_" + std::string( option[ basic::options::OptionKeys::abinitio::log_frags_2 ]() );
	utility::io::ozstream output( filename );
	LoggedFragmentMover& log_frag = dynamic_cast< LoggedFragmentMover& > (*brute_move_large_);
	log_frag.show( output );
	log_frag.clear();
	} // if option[ abinitio::log_frags ].user()
	*/

} // MembraneAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix )

Size
MembraneAbinitio::get_tmh_inserted(core::pose::Pose const & pose) const
{
	return MembraneTopology_from_pose(pose).tmh_inserted();
}
/// @details Pose must already contain a cenlist object or this method will fail.
core::scoring::MembraneTopology const &
MembraneAbinitio::MembraneTopology_from_pose( core::pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;

	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
core::scoring::MembraneTopology &
MembraneAbinitio::nonconst_MembraneTopology_from_pose( core::pose::Pose & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ) {
		return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
	}
	// else
	core::scoring::MembraneTopologyOP membrane_topology( new core::scoring::MembraneTopology );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, membrane_topology );
	return *membrane_topology;
}
/// @details Pose must already contain a cenlist object or this method will fail.
core::scoring::MembraneEmbed const &
MembraneAbinitio::MembraneEmbed_from_pose( core::pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_EMBED;

	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ));
}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
core::scoring::MembraneEmbed &
MembraneAbinitio::nonconst_MembraneEmbed_from_pose( core::pose::Pose & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_EMBED;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ) {
		return *( utility::pointer::static_pointer_cast< core::scoring::MembraneEmbed > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ));
	}
	// else
	core::scoring::MembraneEmbedOP membrane_embed( new core::scoring::MembraneEmbed );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED, membrane_embed );
	return *membrane_embed;
}


} //abinitio
} //protocols
