// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/StepWiseMonteCarlo.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/util.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

//Req'd on WIN32
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>

using namespace core;
using namespace protocols::stepwise::monte_carlo::mover;
using namespace protocols::stepwise::monte_carlo::options;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.StepWiseMonteCarlo" );
using ObjexxFCL::lead_zero_string_of;

//////////////////////////////////////////////////////////////////////////
// StepWiseMonteCarlo -- monte carlo minimization framework for
//  modeler RNA and moves that delete or add residues at chain termini.
//
// 1. This holds the StepWiseMasterMover, and can run the full MonteCarlo or
//     a single move [through the set_move() option, or -move from command line].
// 2. Handles output of 'movie', i.e. monte carlo frames.
// 3. Handles any annealing of scorefxn during monte carlo.
//
//
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace monte_carlo {

//Constructor
StepWiseMonteCarlo::StepWiseMonteCarlo( core::scoring::ScoreFunctionCOP scorefxn_input ):
	scorefxn_input_( scorefxn_input ),
	options_( StepWiseMonteCarloOptionsCOP( new StepWiseMonteCarloOptions ) ), // can be replaced later
	master_mover_( StepWiseMasterMoverOP( new StepWiseMasterMover( scorefxn_input_, options_ ) ) ),
	max_missing_weight_( scorefxn_input->get_weight( scoring::missing_res ) ), // for annealing
	missing_weight_interval_( 0.0 ), // updated below
	missing_weight_( 0.0 ), // can change during run.
	model_tag_( "" ),
	out_path_( "" ),
	out_file_prefix_( "" ),
	movie_file_trial_( "" ),
	movie_file_accepted_( "" ),
	start_time_( clock() )
{
}

//Destructor
StepWiseMonteCarlo::~StepWiseMonteCarlo()
{}

/////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::initialize() {
	initialize_scorefunction();
	master_mover_->set_native_pose( get_native_pose() );
	master_mover_->initialize( scorefxn_, options_ );
	start_time_ = clock();
}

////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::apply( core::pose::Pose & pose ) {
	initialize();

	if ( checkpoint_file_exists() ) {
		pose = pose_from_checkpoint_file();
	}

	if ( pose.size() > 0 && move_.move_type() == NO_MOVE ) show_scores( pose, "Initial score:" );
	if ( master_mover_->do_test_move( move_, pose ) ) {
		show_scores( pose, "After-move score:" );
		return;
	}

	master_mover_->initialize_pose_if_empty( pose );
	if ( !options_->skip_preminimize() ) master_mover_->preminimize_pose( pose );
	do_main_loop( pose );
}

/////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::do_main_loop( pose::Pose & pose ){
	using namespace protocols::moves;
	using namespace core::scoring;

	MonteCarloOP monte_carlo( new MonteCarlo( pose, *scorefxn_, options_->temperature() ) );

	// Oh, if checkpoint_file_exists(), we also want to recover our old low.
	if ( checkpoint_file_exists() ) {
		Pose low = pose_from_checkpoint_file_low();
		monte_carlo->set_lowest_score_pose( low );
	}


	initialize_for_movie( pose );

	Size k( 0 );
	// Restore value of k from pose checkpoint, if possible.
	if ( hasPoseExtraScore( pose, "frame" ) ) {
		k = getPoseExtraScore( pose, "frame" );
	}

	bool success( true );
	Real before_move_score( 0.0 ), after_move_score( 0.0 );
	modeler::switch_focus_among_poses_randomly( pose, scorefxn_ );

	////////////////
	// Main loop
	////////////////

	// If options_->continue_until_none_missing() is true, also require there to be
	// zero missing at end. Caveat: if there are bulge_res, we should account for
	// them appropriately (it's kind of a legacy option, but it still matters)
	Size n_missing = pose.energies().total_energies()[ missing_res ]/ pose.energies().weights()[ missing_res ];
	Size n_bulge = core::pose::full_model_info::const_full_model_info( pose ).rna_bulge_res().size();
	bool none_missing_or_incomplete_okay = ( !options_->continue_until_none_missing() || n_missing == n_bulge );

	while ( k < Size( options_->cycles() ) && none_missing_or_incomplete_okay ) {

		if ( success ) before_move_score = display_progress( pose, k+1 );

		master_mover_->apply( pose );
		success = master_mover_->success();
		if ( !success ) continue;
		k++;

		TR << "After move, modeling: " << get_all_res_list( pose ) << std::endl;
		after_move_score = show_scores( pose, "After-move score:" );

		Real const & proposal_density_ratio = master_mover_->proposal_density_ratio();
		TR << "Score changed from: " << before_move_score << " to  " << after_move_score;
		if ( proposal_density_ratio != 1 ) TR <<  " [proposal_density_ratio: " << proposal_density_ratio << "] ";
		TR << std::endl;
		output_movie( pose, k, "TRIAL", movie_file_trial_ );

		anneal_missing( monte_carlo ); // in use? deprecate?

		monte_carlo->boltzmann( pose, master_mover_->move_type_string(), proposal_density_ratio );

		TR << "Monte Carlo accepted? " << monte_carlo->mc_accepted_string() << std::endl;
		monte_carlo->show_counters();
		output_movie( pose, k, "ACCEPTED", movie_file_accepted_ );

		if ( options_->checkpoint() && k % options_->checkpointing_frequency() == 0 ) {
			remove_checkpoint_file();
			output_checkpoint_file( pose, monte_carlo->lowest_score_pose(), k );
		}
	}
	// Done with this pose, so we can remove the checkpoint file
	if ( options_->checkpoint() )	remove_checkpoint_file();

	if ( options_->recover_low() ) monte_carlo->recover_low( pose );
	show_scores( pose, "Final score:" );

	clearPoseExtraScores( pose );
	if ( options_->save_times() ) setPoseExtraScore( pose, "time", static_cast< Real >( clock() - start_time_ ) / CLOCKS_PER_SEC );
}

// AMW: no need to have k in there because we only want latest
void
StepWiseMonteCarlo::output_checkpoint_file( pose::Pose const & pose, pose::Pose const & lowest, Size const k ) const {
	// I don't think we'll need this. In theory we could store a bunch of checkpoints and restore only the one with biggest frame.
	Pose pose_copy = pose;
	setPoseExtraScore(pose_copy, "frame", k);
	output_to_silent_file( model_tag_ + "_CHECK" /*+ lead_zero_string_of( k, 6 )*/, checkpoint_file_name(), pose_copy, get_native_pose() );

	Pose low_copy = lowest;
	setPoseExtraScore(low_copy, "frame", k);
	output_to_silent_file( model_tag_ + "_CHECK_LOW" /*+ lead_zero_string_of( k, 6 )*/, checkpoint_file_name_low(), low_copy, get_native_pose() );
}

// AMW: no need to have k in there because we only want latest
void
StepWiseMonteCarlo::remove_checkpoint_file() const {
	core::io::silent::remove_silent_file_if_it_exists( checkpoint_file_name() );
}

bool
StepWiseMonteCarlo::checkpoint_file_exists() const {
	return utility::file::file_exists( checkpoint_file_name() );
}

core::pose::Pose
StepWiseMonteCarlo::pose_from_checkpoint_file() const {
	core::pose::Pose pose;
	std::string const checkpoint_file = checkpoint_file_name();
	core::io::silent::SilentFileOptions opts;
	auto sfd = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( opts ) );
	debug_assert( checkpoint_file_exists() );
	sfd->read_file( checkpoint_file );

	core::io::silent::SilentStructOP s( sfd->structure_list()[ 1 ] );
	s->fill_pose( pose );

	return pose;
}

core::pose::Pose
StepWiseMonteCarlo::pose_from_checkpoint_file_low() const {
	core::pose::Pose pose;
	std::string const checkpoint_file = checkpoint_file_name_low();
	core::io::silent::SilentFileOptions opts;
	auto sfd = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( opts ) );
	debug_assert( checkpoint_file_exists() );
	sfd->read_file( checkpoint_file );

	core::io::silent::SilentStructOP s( sfd->structure_list()[ 1 ] );
	s->fill_pose( pose );

	return pose;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::initialize_scorefunction(){
	scorefxn_ = scorefxn_input_->clone();
	scorefxn_->set_weight( scoring::missing_res, 0.0 );
	missing_weight_interval_ = max_missing_weight_ / static_cast<Real>( options_->cycles() );
	missing_weight_ =  missing_weight_interval_;

	if ( options_->rebuild_bulge_mode() ) {
		missing_weight_ = 100.0;
		missing_weight_interval_ = 0.0;
		scorefxn_->set_weight( scoring::missing_res, missing_weight_ );
		scorefxn_->set_weight( scoring::loop_close, 0.0 );
	}

	master_mover_->set_scorefxn( scorefxn_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::initialize_for_movie( pose::Pose const & pose ){
	using namespace utility::file;

	if ( !options_->make_movie() ) return;
	runtime_assert( model_tag_.size() > 0 );
	std::string const dirname =  out_path_ + "movie";
	TR << "Going to make: " << dirname << std::endl;
	create_directory( dirname );

	Pose pose_copy = pose;
	setPoseExtraScore(pose_copy, "frame", 0);

	movie_file_trial_    = out_path_ + "movie/" + model_tag_ + "_trial.out";
	if ( file_exists( movie_file_trial_ ) ) std::remove( movie_file_trial_.c_str() );
	output_movie( pose_copy, 0, "TRIAL", movie_file_trial_ );

	movie_file_accepted_ = out_path_ + "movie/" + model_tag_ + "_accepted.out";
	if ( file_exists( movie_file_accepted_ ) ) std::remove( movie_file_accepted_.c_str() );
	output_movie( pose_copy, 0, "ACCEPTED", movie_file_accepted_ );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::output_movie( pose::Pose const & pose, Size const k, std::string const & tag, std::string const & movie_file ){
	if ( !options_->make_movie() ) return;
	Pose pose_copy = pose;
	setPoseExtraScore(pose_copy, "frame", k);
	output_to_silent_file( tag + "_" + lead_zero_string_of( k, 6 ), movie_file, pose_copy, get_native_pose() );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::anneal_missing( protocols::moves::MonteCarloOP monte_carlo ){
	monte_carlo->change_weight( scoring::missing_res, missing_weight_ );
	missing_weight_ += missing_weight_interval_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseMonteCarlo::display_progress( pose::Pose & pose, Size const cycle_num ){
	using namespace core::pose::full_model_info;
	TR << std::endl << TR.Blue << "Embarking on cycle " << cycle_num << " of " << options_->cycles() << TR.Reset << std::endl;
	Real const before_move_score = show_scores( pose, "Before-move score:" );
	TR << "Modeling: " << get_all_res_list( pose ) << std::endl;
	if ( missing_weight_ != 0.0 ) TR << "Weight of missing residue score term is " << missing_weight_ << std::endl;
	return before_move_score;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseMonteCarlo::show_scores( core::pose::Pose & pose,
	std::string const & tag ){
	if ( options_->verbose_scores() ) {
		TR << tag << " " << ( *scorefxn_ )( pose ) << std::endl;
		scorefxn_->show( TR, pose );
	}
	return (*scorefxn_)( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::set_options( options::StepWiseMonteCarloOptionsCOP options ){
	options_ = options;
	master_mover_->set_options( options );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::set_submotif_library( monte_carlo::submotif::SubMotifLibraryCOP setting ) {
	master_mover_->set_submotif_library( setting );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::set_native_pose( PoseCOP pose ){
	Mover::set_native_pose( pose );
	master_mover_->set_native_pose( pose );
}

/// @brief Helper that contains the formula for naming checkpoint files.
std::string
StepWiseMonteCarlo::checkpoint_file_name() const {
	return out_path_ + out_file_prefix_ + ".checkpoint"; // do not end in .out, or will get collated by easy_cat.py
}

/// @brief Helper that contains the formula for naming checkpoint files.
std::string
StepWiseMonteCarlo::checkpoint_file_name_low() const {
	return out_path_ + out_file_prefix_ + ".checkpoint_low"; // do not end in .out, or will get collated by easy_cat.py
}

} //monte_carlo
} //stepwise
} //protocols
