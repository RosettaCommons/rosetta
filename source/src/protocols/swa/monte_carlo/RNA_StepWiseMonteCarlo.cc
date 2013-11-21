// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/RNA_StepWiseMonteCarlo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/swa/monte_carlo/RNA_StepWiseMonteCarlo.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/swa/monte_carlo/SWA_MonteCarloUtil.hh>
#include <protocols/swa/monte_carlo/RNA_AddMover.hh>
#include <protocols/swa/monte_carlo/RNA_DeleteMover.hh>
#include <protocols/swa/monte_carlo/RNA_AddOrDeleteMover.hh>
#include <protocols/swa/monte_carlo/RNA_ResampleMover.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>

using namespace core;

static numeric::random::RandomGenerator RG(2391021);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.swa.monte_carlo.RNA_StepWiseMonteCarlo" );
using ObjexxFCL::string_of;

//////////////////////////////////////////////////////////////////////////
// StepWiseMonteCarlo -- monte carlo minimization framework for
//  sampling RNA and moves that delete or add residues at chain termini.
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace swa {
namespace monte_carlo {

//Constructor
	RNA_StepWiseMonteCarlo::RNA_StepWiseMonteCarlo( core::scoring::ScoreFunctionOP scorefxn, core::pose::PoseOP native_pose, core::Real constraint_x0, core::Real constraint_tol ):
	scorefxn_( scorefxn ),
	verbose_scores_( false ),
	use_phenix_geo_( true ),
	skip_deletions_( false ),
	erraser_( true ),
	allow_internal_moves_( false ),
	num_random_samples_( 20 ),
	cycles_( 500 ),
	add_delete_frequency_( 0.5 ),
	minimize_single_res_frequency_( 0.0 ),
	switch_focus_frequency_( 0.5 ),
	just_min_after_mutation_frequency_( 0.5 ),
	temperature_( 1.0 ),
	native_pose_( native_pose ),
	constraint_x0_( constraint_x0 ),
	constraint_tol_( constraint_tol )
{
	using namespace core::scoring;
	//rmsd_weight_ = scorefxn_->get_weight( coordinate_constraint );
	max_missing_weight_ = scorefxn_->get_weight( missing_res );
	chainbreak_weight_ = scorefxn_->get_weight( linear_chainbreak );
}


//Destructor
RNA_StepWiseMonteCarlo::~RNA_StepWiseMonteCarlo()
{}

////////////////////////////////////////////////////////////////////////////////
void
RNA_StepWiseMonteCarlo::apply( core::pose::Pose & pose ) {

	using namespace protocols::moves;
	using namespace protocols::swa;
	using namespace core::scoring;

	initialize_movers();
	//scorefxn_->set_weight( coordinate_constraint, rmsd_weight_ );
	scorefxn_->set_weight( missing_res, 0.0 );
	ScoreFunctionOP temp_score = scorefxn_->clone();
	MonteCarloOP monte_carlo_ = new MonteCarlo( pose, *temp_score, temperature_ );
	//scorefxn_->set_weight( coordinate_constraint, 0.0 );
	show_scores( pose, "Initial score:" );

	Size k( 1 );
	std::string move_type;
	bool success( true );
	Real const missing_weight_interval = max_missing_weight_ / cycles_;
	//missing_weight_interval /= cycles_;
	Real missing_weight = missing_weight_interval;

	while (  k <= cycles_ ){
		//scorefxn_->set_weight( missing_res, missing_weight );
		//clear_constraints_recursively( pose );

		if ( success ) {
			TR << std::endl << TR.Blue << "Embarking on cycle " << k << " of " << cycles_ << TR.Reset << std::endl;
			show_scores( pose, "Before-move score:" );
			TR << "Weight of missing residue score term is " << missing_weight << std::endl;
		}

		if ( RG.uniform() < switch_focus_frequency_ ) switch_focus_among_poses_randomly( pose );

		bool const minimize_single_res = ( RG.uniform() <= minimize_single_res_frequency_ );

		if ( RG.uniform() < add_delete_frequency_ ) {
				rna_add_or_delete_mover_->set_minimize_single_res( minimize_single_res );
				success = rna_add_or_delete_mover_->apply( pose, move_type );
		} else {
			// later make this an actual class!
				rna_resample_mover_->set_minimize_single_res( minimize_single_res );
				success = rna_resample_mover_->apply( pose, move_type );
		}

		if ( !success ) continue;
		k++;
		show_scores( pose, "After-move score:" );
		
		scoring::EMapVector & energy_map = pose.energies().total_energies();
		Real const chainbreak_score = energy_map[linear_chainbreak];
		
		if ( chainbreak_score > 0.01 ) {
			monte_carlo_->change_weight( linear_chainbreak, 100 );
		}
		
		monte_carlo_->change_weight( missing_res, missing_weight );
		missing_weight += missing_weight_interval;
		if ( minimize_single_res ) move_type += "-minsngl";
		
//		if ( native_pose_ ) {
//			superimpose_recursively_and_add_constraints( pose, *native_pose_ );
//		}
		
		monte_carlo_->boltzmann( pose, move_type );
		
		monte_carlo_->change_weight( linear_chainbreak, chainbreak_weight_);

		// following can be removed later.
		TR << "Monte Carlo accepted? " << monte_carlo_->mc_accepted_string() << std::endl;
		monte_carlo_->show_counters();

	}

	monte_carlo_->recover_low( pose );
	scorefxn_->set_weight( missing_res, max_missing_weight_ );
	//scorefxn_->set_weight( coordinate_constraint, rmsd_weight_ );
	show_scores( pose, "Final score:" );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StepWiseMonteCarlo::initialize_movers(){

	using namespace protocols::swa::rna;

	// used in all movers.
	StepWiseRNA_ModelerOP stepwise_rna_modeler = new StepWiseRNA_Modeler( scorefxn_ );
	stepwise_rna_modeler->set_choose_random( true );
	stepwise_rna_modeler->set_force_centroid_interaction( true );
	stepwise_rna_modeler->set_use_phenix_geo( use_phenix_geo_ );
	stepwise_rna_modeler->set_kic_sampling_if_relevant( erraser_ );
	stepwise_rna_modeler->set_num_random_samples( num_random_samples_ );
	stepwise_rna_modeler->set_num_pose_minimize( 1 );

	// maybe RNA_AddMover could just hold a copy of RNA_ResampleMover...
	rna_add_mover_ = new RNA_AddMover( scorefxn_, native_pose_, constraint_x0_, constraint_tol_ );
	rna_add_mover_->set_start_added_residue_in_aform( false );
	rna_add_mover_->set_presample_added_residue(  true );
	rna_add_mover_->set_presample_by_swa(  true );
	rna_add_mover_->set_stepwise_rna_modeler( stepwise_rna_modeler->clone_modeler() );

	rna_delete_mover_ = new RNA_DeleteMover( native_pose_, constraint_x0_, constraint_tol_ );
	// following will be used to minimize after deletion.
	rna_delete_mover_->set_stepwise_rna_modeler( stepwise_rna_modeler->clone_modeler() );

	rna_add_or_delete_mover_ = new RNA_AddOrDeleteMover( rna_add_mover_, rna_delete_mover_ );
	rna_add_or_delete_mover_->set_sample_res( sample_res_ );
	rna_add_or_delete_mover_->set_skip_deletions( skip_deletions_ ); // for testing only

	rna_resample_mover_ = new RNA_ResampleMover( stepwise_rna_modeler->clone_modeler(), native_pose_, constraint_x0_, constraint_tol_ );
	rna_resample_mover_->set_just_min_after_mutation_frequency( just_min_after_mutation_frequency_ );
	rna_resample_mover_->set_allow_internal_moves( allow_internal_moves_ );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_StepWiseMonteCarlo::switch_focus_among_poses_randomly( pose::Pose & pose ) const {

	using namespace core::pose;
	using namespace core::pose::full_model_info;
	using namespace protocols::swa;

	Size const num_other_poses = const_full_model_info( pose ).other_pose_list().size();
	Size const focus_pose_idx = RG.random_range( 0, num_other_poses );
	if ( focus_pose_idx == 0 ) return false;

	Real const score_before_switch_focus = (*scorefxn_)( pose );
	TR.Debug << "SWITCHING FOCUS! SWITCHING FOCUS! SWITCHING FOCUS! SWITCHING FOCUS! to: " << focus_pose_idx << std::endl;
	switch_focus_to_other_pose( pose, focus_pose_idx );
	Real const score_after_switch_focus = (*scorefxn_)( pose );

	// originally set threshold at 0.001, but triggered rare errors. At some point worth tracking down...
	if (  std::abs( score_before_switch_focus - score_after_switch_focus ) > 0.10 ){
		utility_exit_with_message( "Energy change after switching pose focus: " + string_of( score_before_switch_focus ) + " to " +string_of( score_after_switch_focus ) );
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_StepWiseMonteCarlo::show_scores( core::pose::Pose & pose,
																		 std::string const tag ){

	TR << tag << " " << ( *scorefxn_ )( pose ) << std::endl;
	if ( verbose_scores_ ) {
		scorefxn_->show( TR, pose );
	}
}


} //monte_carlo
} //swa
} //protocols
