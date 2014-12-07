// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/StepWiseMonteCarlo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMover.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.hh>
#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.hh>
#include <protocols/stepwise/monte_carlo/mover/ResampleMover.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/precomputed/PrecomputedLibraryMover.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
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

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

//Req'd on WIN32
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>

using namespace core;
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::monte_carlo::mover;
using namespace protocols::stepwise::monte_carlo::rna;

static thread_local basic::Tracer TR( "protocols.stepwise.monte_carlo.StepWiseMonteCarlo" );
using ObjexxFCL::lead_zero_string_of;

//////////////////////////////////////////////////////////////////////////
// StepWiseMonteCarlo -- monte carlo minimization framework for
//  modeler RNA and moves that delete or add residues at chain termini.
//
// This is the master Mover, and can run the full MonteCarlo or
//  a single move [through the set_move() option].
//
// This also handles preminimizing. But it may be better to move the
//  preminimizer out of here -- it can have a different setting
//  (e.g. use_packer_instead_of_rotamer_trials) and perhaps should
//  not output to silent file if there's going to be a subsequent move.
//
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace monte_carlo {

//Constructor
StepWiseMonteCarlo::StepWiseMonteCarlo( core::scoring::ScoreFunctionCOP scorefxn_input ):
	scorefxn_input_( scorefxn_input ),
	options_( options::StepWiseMonteCarloOptionsCOP( options::StepWiseMonteCarloOptionsOP( new options::StepWiseMonteCarloOptions ) ) ), // can be replaced later
	minimize_single_res_( false ), // changes during run
	max_missing_weight_( scorefxn_input->get_weight( scoring::missing_res ) ), // for annealing
	missing_weight_interval_( 0.0 ), // updated below
	missing_weight_( 0.0 ), // can change during run.
	model_tag_( "" ),
	out_path_( "" ),
	movie_file_trial_( "" ),
	movie_file_accepted_( "" ),
	enumerate_( false ),
	do_preminimize_move_( false ),
	start_time_( clock() )
{
}

//Destructor
StepWiseMonteCarlo::~StepWiseMonteCarlo()
{}

////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::apply( core::pose::Pose & pose ) {
	initialize();

	if ( pose.total_residue() > 0 ) show_scores( pose, "Initial score:" );
	if ( do_test_move( pose ) ) return;

	initialize_pose_if_empty( pose );
	preminimize_pose( pose );
	do_main_loop( pose );
}

/////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::initialize() {
	initialize_scorefunction();
	initialize_movers();
	start_time_ = clock();
}

/////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::do_main_loop( pose::Pose & pose ){
	using namespace protocols::moves;
	using namespace core::scoring;

	MonteCarloOP monte_carlo( new MonteCarlo( pose, *scorefxn_, options_->temperature() ) );
	initialize_for_movie( pose );

	Size k( 0 );
	std::string move_type;
	bool success( true );
	Real before_move_score( 0.0 ), after_move_score( 0.0 );
	switch_focus_among_poses_randomly( pose, scorefxn_ );

	////////////////
	// Main loop
	////////////////
	while ( k < Size( options_->cycles() ) ){

		if ( success ) before_move_score = display_progress( pose, k+1 );

		if ( numeric::random::rg().uniform() < options_->switch_focus_frequency() ) switch_focus_among_poses_randomly( pose );

		set_minimize_single_res( numeric::random::rg().uniform() <= options_->minimize_single_res_frequency() );
		if ( numeric::random::rg().uniform() < options_->add_delete_frequency() ) {
			success = add_or_delete_mover_->apply( pose, move_type );
		} else {
			success = resample_mover_->apply( pose, move_type );
		}
		if ( minimize_single_res_ ) move_type += "-minsngl";

		if ( !success ) continue;
		k++;

		TR << "After move, modeling: " << get_all_res_list( pose ) << std::endl;
		after_move_score = show_scores( pose, "After-move score:" );
		TR << "Score changed from: " << before_move_score << " to  " << after_move_score << std::endl;
		output_movie( pose, k, "TRIAL", movie_file_trial_ );

		anneal_missing( monte_carlo );
		monte_carlo->boltzmann( pose, move_type );

		TR << "Monte Carlo accepted? " << monte_carlo->mc_accepted_string() << std::endl;
		monte_carlo->show_counters();
		output_movie( pose, k, "ACCEPTED", movie_file_accepted_ );
	}

	if ( options_->recover_low() ) monte_carlo->recover_low( pose );
	show_scores( pose, "Final score:" );

	clearPoseExtraScores( pose );
	if ( options_->save_times() ) setPoseExtraScore( pose, "time", static_cast< Real >( clock() - start_time_ ) / CLOCKS_PER_SEC );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::initialize_scorefunction(){
	scorefxn_ = scorefxn_input_->clone();
	scorefxn_->set_weight( scoring::missing_res, 0.0 );
	missing_weight_interval_ = max_missing_weight_ / static_cast<Real>( options_->cycles() );
	missing_weight_ =  missing_weight_interval_;

	if ( options_->rebuild_bulge_mode() ){
		missing_weight_ = 100.0;
		missing_weight_interval_ = 0.0;
		scorefxn_->set_weight( scoring::missing_res, missing_weight_ );
		scorefxn_->set_weight( scoring::loop_close, 0.0 );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::initialize_movers(){

	using namespace protocols::stepwise::modeler::rna;

	stepwise_modeler_ = setup_unified_stepwise_modeler();
	stepwise_modeler_->set_native_pose( get_native_pose() );

	// maybe AddMover could just hold a copy of ResampleMover...
	add_mover_ = AddMoverOP( new AddMover( scorefxn_ ) );
	add_mover_->set_native_pose( get_native_pose() );
	add_mover_->set_start_added_residue_in_aform( false );
	add_mover_->set_presample_added_residue(  true );
	add_mover_->set_presample_by_swa(  true );
	add_mover_->set_stepwise_modeler( stepwise_modeler_->clone_modeler() );

	delete_mover_ = DeleteMoverOP( new DeleteMover );
	delete_mover_->set_native_pose( get_native_pose() );
	delete_mover_->set_stepwise_modeler( stepwise_modeler_->clone_modeler() );
	delete_mover_->set_options( options_ );

	from_scratch_mover_ = FromScratchMoverOP( new FromScratchMover );
	from_scratch_mover_->set_native_pose( get_native_pose() );
	from_scratch_mover_->set_stepwise_modeler( stepwise_modeler_->clone_modeler() );

	add_or_delete_mover_ = AddOrDeleteMoverOP( new AddOrDeleteMover( add_mover_, delete_mover_, from_scratch_mover_ ) );
	add_or_delete_mover_->set_options( options_ );

	resample_mover_ = ResampleMoverOP( new ResampleMover( stepwise_modeler_->clone_modeler() ) );
	resample_mover_->set_native_pose( get_native_pose() );
	resample_mover_->set_options( options_ );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// used in all movers.
modeler::StepWiseModelerOP
StepWiseMonteCarlo::setup_unified_stepwise_modeler(){
 	using namespace modeler;
	using namespace modeler::rna;
 	using namespace modeler::protein;
 	using namespace modeler::precomputed;

 	StepWiseModelerOP stepwise_modeler( new StepWiseModeler( scorefxn_ ) );
	protocols::stepwise::modeler::options::StepWiseModelerOptionsOP options = options_->setup_modeler_options();
	if ( enumerate_ ) options->set_choose_random( false );
	if ( do_preminimize_move_ ) options->set_use_packer_instead_of_rotamer_trials( true ); // for proteins.
	stepwise_modeler->set_options( options );

	if ( options_->use_precomputed_library() ){
		PrecomputedLibraryMoverOP precomputed_library_mover( new PrecomputedLibraryMover );
		stepwise_modeler->set_precomputed_library_mover( precomputed_library_mover );
	}
 	return stepwise_modeler;
 }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::initialize_pose_if_empty( pose::Pose & pose ){
	if ( pose.total_residue() > 0 ) return;
	runtime_assert( options_->from_scratch_frequency() > 0.0 );
	add_or_delete_mover_->apply( pose );
	runtime_assert( pose.total_residue() > 0 );
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

	movie_file_trial_    = out_path_ + "movie/" + model_tag_ + "_trial.out";
	if ( file_exists( movie_file_trial_ ) ) std::remove( movie_file_trial_.c_str() );
	output_movie( pose, 0, "TRIAL", movie_file_trial_ );

	movie_file_accepted_ = out_path_ + "movie/" + model_tag_ + "_accepted.out";
	if ( file_exists( movie_file_accepted_ ) ) std::remove( movie_file_accepted_.c_str() );
	output_movie( pose, 0, "ACCEPTED", movie_file_accepted_ );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::output_movie( pose::Pose const & pose, Size const k, std::string const tag, std::string const & movie_file ){
	if ( !options_->make_movie() ) return;
	Pose pose_copy = pose;
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

//////////////////////////////////////////////////////////////////////////////////////
std::string
StepWiseMonteCarlo::get_all_res_list( pose::Pose & pose ){
	std::string out_string;
	out_string = make_tag_with_dashes( get_res_list_from_full_model_info_const( pose ) );
	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	if ( other_pose_list.size() == 0 ) return out_string;
	out_string += " [";
	for ( Size n = 1; n <= other_pose_list.size(); n++ ){
		out_string += get_all_res_list( *other_pose_list[n] );
		if ( n < other_pose_list.size() ) out_string += "; ";
	}
	out_string += "]";
	return out_string;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseMonteCarlo::show_scores( core::pose::Pose & pose,
																		 std::string const tag ){
	if ( options_->verbose_scores() ) {
		TR << tag << " " << ( *scorefxn_ )( pose ) << std::endl;
		scorefxn_->show( TR, pose );
	}
	return (*scorefxn_)( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::set_minimize_single_res( bool const minimize_single_res ){
	add_or_delete_mover_->set_minimize_single_res( minimize_single_res );
	resample_mover_->set_minimize_single_res( minimize_single_res );
	minimize_single_res_ = minimize_single_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::set_options( options::StepWiseMonteCarloOptionsCOP options ){
	options_ = options;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseMonteCarlo::do_test_move( pose::Pose & pose ){
	if ( do_preminimize_move_ ) preminimize_pose( pose );
	if ( move_.move_type() == NO_MOVE ) return do_preminimize_move_;
	if ( move_.move_type() == ADD || move_.move_type() == DELETE || move_.move_type() == FROM_SCRATCH ){
		add_or_delete_mover_->apply( pose, move_ );
	} else {
		resample_mover_->apply( pose, move_ );
	}
	show_scores( pose, "After-move score:" );
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMonteCarlo::preminimize_pose( pose::Pose & pose ) {
	fix_up_residue_type_variants( pose );
	stepwise_modeler_->set_moving_res_and_reset( 0 );
	stepwise_modeler_->set_working_prepack_res( get_all_residues( pose ) );
	stepwise_modeler_->set_working_minimize_res( get_moving_res_from_full_model_info( pose ) );
	stepwise_modeler_->apply( pose );
	utility::vector1< PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
	for ( Size n = 1; n <= other_pose_list.size(); n++ ){
		preminimize_pose( *( other_pose_list[ n ] ) );
	}
}


/////////////////////////////////////////////////////////
AddOrDeleteMoverOP
StepWiseMonteCarlo::add_or_delete_mover(){ return add_or_delete_mover_; }

/////////////////////////////////////////////////////////
// Called by build_full_model() in stepwise/monte_carlo/util.cc
void
StepWiseMonteCarlo::build_full_model( pose::Pose const & start_pose, Pose & full_model_pose ){
	using namespace options;
	full_model_pose = start_pose;
	runtime_assert( options_->skip_deletions() ); // totally inelegant, must be set outside.
	initialize();
	add_or_delete_mover_->set_choose_random( false );
	add_mover_->set_start_added_residue_in_aform( true  );
	add_mover_->set_presample_added_residue(      false );

	std::string move_type;
	while( add_or_delete_mover_->apply( full_model_pose, move_type ) ){
		TR.Debug << "Building full model: " << move_type << std::endl;
	}
}



} //monte_carlo
} //stepwise
} //protocols
