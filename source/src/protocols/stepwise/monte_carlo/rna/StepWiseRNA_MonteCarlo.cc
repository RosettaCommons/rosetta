// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarlo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloUtil.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_AddMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_DeleteMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_FromScratchMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_AddOrDeleteMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_ResampleMover.hh>
#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Modeler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

using namespace core;

static numeric::random::RandomGenerator RG(2391021);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.stepwise.monte_carlo.StepWiseRNA_MonteCarlo" );
using ObjexxFCL::lead_zero_string_of;

//////////////////////////////////////////////////////////////////////////
// StepWiseMonteCarlo -- monte carlo minimization framework for
//  sampling RNA and moves that delete or add residues at chain termini.
//
// This is the master Mover.
//
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

//Constructor
StepWiseRNA_MonteCarlo::StepWiseRNA_MonteCarlo( core::scoring::ScoreFunctionCOP scorefxn_input ):
	scorefxn_input_( scorefxn_input ),
	options_( new StepWiseRNA_MonteCarloOptions ), // can be replaced later
	minimize_single_res_( false ), // changes during run
	max_missing_weight_( scorefxn_input->get_weight( scoring::missing_res ) ), // for annealing
	missing_weight_interval_( 0.0 ), // updated below
	missing_weight_( 0.0 ), // can change during run.
	model_tag_( "" ),
	out_path_( "" ),
	movie_file_trial_( "" ),
	movie_file_accepted_( "" ),
	enumerate_( false )
{
}

//Destructor
StepWiseRNA_MonteCarlo::~StepWiseRNA_MonteCarlo()
{}

////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::apply( core::pose::Pose & pose ) {

	initialize_scorefunction();
	initialize_movers();
	initialize_pose_if_empty( pose );

	if ( do_test_move( pose ) ) return;
	do_main_loop( pose );
}

/////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::do_main_loop( pose::Pose & pose ){
	using namespace protocols::moves;
	using namespace core::scoring;

	MonteCarloOP monte_carlo = new MonteCarlo( pose, *scorefxn_, options_->temperature() );
	show_scores( pose, "Initial score:" );
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

		if ( RG.uniform() < options_->switch_focus_frequency() ) switch_focus_among_poses_randomly( pose );

		set_minimize_single_res( RG.uniform() <= options_->minimize_single_res_frequency() );
		if ( RG.uniform() < options_->add_delete_frequency() ) {
			success = rna_add_or_delete_mover_->apply( pose, move_type );
		} else {
			success = rna_resample_mover_->apply( pose, move_type );
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

	monte_carlo->recover_low( pose );
	show_scores( pose, "Final score:" );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::initialize_scorefunction(){
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
StepWiseRNA_MonteCarlo::initialize_movers(){

	using namespace protocols::stepwise::sampling::rna;

	// used in all movers.
	StepWiseRNA_ModelerOP stepwise_rna_modeler = new StepWiseRNA_Modeler( scorefxn_ );
	StepWiseRNA_ModelerOptionsOP modeler_options = options_->setup_modeler_options();
	stepwise_rna_modeler->set_options( modeler_options );
	stepwise_rna_modeler->set_minimizer_extra_minimize_res( options_->extra_minimize_res() );
	stepwise_rna_modeler->set_syn_chi_res_list( options_->syn_chi_res_list() );
	stepwise_rna_modeler->set_terminal_res( options_->terminal_res() );
	if ( enumerate_ ) { // SWA-like.
 		modeler_options->set_choose_random( false );
		modeler_options->set_num_pose_minimize( 108 );
	}

	// maybe RNA_AddMover could just hold a copy of RNA_ResampleMover...
	rna_add_mover_ = new RNA_AddMover( scorefxn_ );
	rna_add_mover_->set_native_pose( get_native_pose() );
	rna_add_mover_->set_start_added_residue_in_aform( false );
	rna_add_mover_->set_presample_added_residue(  true );
	rna_add_mover_->set_presample_by_swa(  true );
	rna_add_mover_->set_stepwise_rna_modeler( stepwise_rna_modeler->clone_modeler() );
	rna_add_mover_->set_constraint_x0( options_->constraint_x0() );
	rna_add_mover_->set_constraint_tol( options_->constraint_tol() );

	rna_delete_mover_ = new RNA_DeleteMover;
	rna_delete_mover_->set_native_pose( get_native_pose() );
	rna_delete_mover_->set_stepwise_rna_modeler( stepwise_rna_modeler->clone_modeler() );
	rna_delete_mover_->set_options( options_ );

	rna_from_scratch_mover_ = new RNA_FromScratchMover;
	rna_from_scratch_mover_->set_native_pose( get_native_pose() );
	rna_from_scratch_mover_->set_stepwise_rna_modeler( stepwise_rna_modeler->clone_modeler() );

	rna_add_or_delete_mover_ = new RNA_AddOrDeleteMover( rna_add_mover_, rna_delete_mover_, rna_from_scratch_mover_ );
	rna_add_or_delete_mover_->set_options( options_ );

	rna_resample_mover_ = new RNA_ResampleMover( stepwise_rna_modeler->clone_modeler() );
	rna_resample_mover_->set_native_pose( get_native_pose() );
	rna_resample_mover_->set_options( options_ );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::initialize_pose_if_empty( pose::Pose & pose ){
	if ( pose.total_residue() > 0 ) return;
	runtime_assert( options_->from_scratch_frequency() > 0.0 );
	rna_add_or_delete_mover_->apply( pose );
	runtime_assert( pose.total_residue() > 0 );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::initialize_for_movie( pose::Pose const & pose ){
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
StepWiseRNA_MonteCarlo::output_movie( pose::Pose const & pose, Size const k, std::string const tag, std::string const & movie_file ){
	if ( !options_->make_movie() ) return;
	Pose pose_copy = pose;
	output_to_silent_file( tag + "_" + lead_zero_string_of( k, 6 ), movie_file, pose_copy, get_native_pose() );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::anneal_missing( protocols::moves::MonteCarloOP monte_carlo ){
	monte_carlo->change_weight( scoring::missing_res, missing_weight_ );
	missing_weight_ += missing_weight_interval_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseRNA_MonteCarlo::display_progress( pose::Pose & pose, Size const cycle_num ){
	using namespace core::pose::full_model_info;
	TR << std::endl << TR.Blue << "Embarking on cycle " << cycle_num << " of " << options_->cycles() << TR.Reset << std::endl;
	Real const before_move_score = show_scores( pose, "Before-move score:" );
	TR << "Modeling: " << get_all_res_list( pose ) << std::endl;
	if ( missing_weight_ != 0.0 ) TR << "Weight of missing residue score term is " << missing_weight_ << std::endl;
	return before_move_score;
}

//////////////////////////////////////////////////////////////////////////////////////
std::string
StepWiseRNA_MonteCarlo::get_all_res_list( pose::Pose & pose ){
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
StepWiseRNA_MonteCarlo::show_scores( core::pose::Pose & pose,
																		 std::string const tag ){
	if ( options_->verbose_scores() ) {
		TR << tag << " " << ( *scorefxn_ )( pose ) << std::endl;
		scorefxn_->show( TR, pose );
	}
	return (*scorefxn_)( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::set_minimize_single_res( bool const minimize_single_res ){
	rna_add_or_delete_mover_->set_minimize_single_res( minimize_single_res );
	rna_resample_mover_->set_minimize_single_res( minimize_single_res );
	minimize_single_res_ = minimize_single_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_MonteCarlo::set_options( StepWiseRNA_MonteCarloOptionsCOP options ){
	options_ = options;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_MonteCarlo::do_test_move( pose::Pose & pose ){
	if ( move_.move_type() == NO_MOVE ) return false;
	if ( move_.move_type() == ADD || move_.move_type() == DELETE ){
		rna_add_or_delete_mover_->apply( pose, move_ );
		return true;
	} else {
		rna_resample_mover_->apply( pose, move_ );
		return true;
	}
}


} //rna
} //monte_carlo
} //stepwise
} //protocols
