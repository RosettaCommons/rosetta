// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/RNA_FragmentMonteCarlo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/RNA_FragmentMonteCarlo.hh>
#include <protocols/farna/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/farna/RNA_FragmentMover.hh>
#include <protocols/farna/RNA_Fragments.hh>
#include <protocols/farna/RNA_ChunkLibrary.hh>
#include <protocols/farna/RNA_JumpLibrary.hh>
#include <protocols/farna/RNA_LoopCloser.hh>
#include <protocols/farna/RNA_LibraryManager.hh>
#include <protocols/farna/RNA_Minimizer.hh>
#include <protocols/farna/RNA_Relaxer.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/RNA_SecStructInfo.hh>
#include <protocols/farna/FullAtomRNA_Fragments.hh>
#include <protocols/farna/BasePairStepLibrary.hh>
#include <protocols/farna/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/copydofs/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <numeric/random/random.hh>
#include <basic/database/open.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.RNA_FragmentMonteCarlo" );

using namespace core;
using namespace ObjexxFCL::format; // AUTO USING NS
using core::io::pdb::dump_pdb;

namespace protocols {
namespace farna {

//Constructor
RNA_FragmentMonteCarlo::RNA_FragmentMonteCarlo( RNA_FragmentMonteCarloOptionsCOP options ):
	options_( options ),
	monte_carlo_cycles_max_default_( 100000 ),
	monte_carlo_cycles_( 0 ), // will be updated later based on options
	rounds_( 0 ), // will be updated later based on options
	frag_size_( 3 ),
	do_close_loops_( true ), // will be updated later based on options
	refine_pose_( false ),
	jump_change_frequency_( 0.1 ), //  maybe updated based on options, or if rigid-body sampling
	lores_score_early_( 0.0 ),
	lores_score_final_( 0.0 ),
	chunk_coverage_( 0.0 ) // will be updated later
{}

//Destructor
RNA_FragmentMonteCarlo::~RNA_FragmentMonteCarlo()
{}

////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::initialize( pose::Pose & pose ) {
	initialize_libraries( pose );
	initialize_movers();
	initialize_score_functions();
	initialize_parameters();
}

////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::initialize_parameters() {
	// parameters for run.
	jump_change_frequency_ = options_->jump_change_frequency(); // might get overwritten if rigid_body movement, tested later.
	rounds_ = refine_pose_ ? 1 : options_->rounds();
}

////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// This sets up StructureParameters, AllowInsert, and various libraries (fragments, 'chunks', etc.)
//
// It is extremely sensitive to order of operations, unfortunately.
// There are two modes:
//
//  1. De novo (standard FARFAR). This function then has to do some sensitive work
//      on the input pose, which is assumed to be unfolded, and have no virtual phosphates.
//
//  2. Refine_pose mode -- accepts the pose as is. In this mode, the pose should actually be
//      const -- need some way to check this.
//
///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::initialize_libraries( pose::Pose & pose ) {

	rna_structure_parameters_ = RNA_StructureParametersOP( new RNA_StructureParameters );

	// all jumping, secondary structure, base pair constraint, allow_insert information
	// will be stored in a .prm file.
	rna_structure_parameters_->set_bps_moves( options_->bps_moves() );
	if ( refine_pose_ ) {
		if ( get_rna_secstruct( pose ).size() != pose.total_residue() ) set_rna_secstruct( pose, std::string( pose.total_residue(), 'X' ) );
		rna_structure_parameters_->initialize_from_pose( pose );
		rna_structure_parameters_->set_jump_library( RNA_LibraryManager::get_instance()->rna_jump_library_cop( options_->jump_library_file() ) );
	} else {
		// assumes folding from scratch:
		//  -- pose is assumed to have no variants (!), including any virtual phosphates or chainbreaks.
		//  -- following actually adds virtual anchor to pose, updates pose secstruct
		//  -- also: handles read in of params file.
		rna_structure_parameters_->initialize_for_de_novo_protocol( pose, options_->rna_params_file(), options_->jump_library_file(), options_->ignore_secstruct() );
	}
	rna_structure_parameters_->set_root_at_first_rigid_body( options_->root_at_first_rigid_body() );
	rna_structure_parameters_->set_suppress_bp_constraint( options_->suppress_bp_constraint() );

	rna_chunk_library_ = protocols::farna::RNA_ChunkLibraryOP( new RNA_ChunkLibrary( rna_structure_parameters_->allow_insert(),
		options_->chunk_pdb_files(), options_->chunk_silent_files(),
		pose, options_->input_res() ) );

	if ( options_->bps_moves() ) {
		rna_chunk_library_->setup_base_pair_step_chunks( pose, rna_structure_parameters_->get_canonical_base_pair_steps(),
			RNA_LibraryManager::get_instance()->canonical_base_pair_step_library() );
		rna_chunk_library_->setup_base_pair_step_chunks( pose, rna_structure_parameters_->get_noncanonical_base_pair_steps(),
			RNA_LibraryManager::get_instance()->general_base_pair_step_library() );
	}

	// following is used to figure out default frequency for chunk insertions.
	chunk_coverage_ = rna_chunk_library_->chunk_coverage();


	rna_structure_parameters_->allow_insert()->and_allow_insert( rna_chunk_library_->allow_insert() );

	if ( !refine_pose_ ) {
		// if input pose is setup from scratch, this is the place where we add variants, and update allow_insert.
		// this still does not include chainbreak variants which are set up (randomly) in the main loop.
		// there's potentially a much better way to handle all this -- setup a pose outside, and then skip this.
		std::cout << pose.annotated_sequence() << std::endl;

		// note this crazy order.
		rna_structure_parameters_->setup_virtual_phosphate_variants( pose ); // needed to refreeze virtual phosphates!
	}
	rna_chunk_library_->set_allow_insert( rna_structure_parameters_->allow_insert() );

}

////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::initialize_movers() {

	protocols::farna::RNA_Fragments const & all_rna_fragments_( RNA_LibraryManager::get_instance()->rna_fragment_library( options_->all_rna_fragments_file() ) );
	rna_fragment_mover_ = protocols::farna::RNA_FragmentMoverOP( new RNA_FragmentMover( all_rna_fragments_, rna_structure_parameters_->allow_insert() ) );

	rna_loop_closer_ = protocols::farna::RNA_LoopCloserOP( new protocols::farna::RNA_LoopCloser );

	rna_minimizer_ = protocols::farna::RNA_MinimizerOP( new RNA_Minimizer );
	rna_minimizer_->set_score_function( hires_scorefxn_ );
	rna_minimizer_->set_allow_insert( rna_structure_parameters_->allow_insert() );
	rna_minimizer_->vary_bond_geometry( options_->vary_bond_geometry() );
	rna_minimizer_->set_extra_minimize_res( options_->extra_minimize_res() );
	rna_minimizer_->set_extra_minimize_chi_res( options_->extra_minimize_chi_res() );
	rna_minimizer_->set_move_first_rigid_body( options_->move_first_rigid_body() );
	rna_minimizer_->use_coordinate_constraints( options_->minimizer_use_coordinate_constraints() );

	rna_relaxer_ = protocols::farna::RNA_RelaxerOP( new RNA_Relaxer( rna_fragment_mover_, rna_minimizer_) );
	rna_relaxer_->simple_rmsd_cutoff_relax( options_->simple_rmsd_cutoff_relax() );
}

////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::initialize_score_functions() {
	// scorefxns
	using namespace core::scoring;
	working_denovo_scorefxn_ = denovo_scorefxn_->clone();
	// RNA high-resolution score function + rna_chem_shift term
	if ( options_->use_chem_shift_data() ) {
		Real const CS_weight = 4.0; //hard-coded to 4.0 based on CS-ROSETTA-RNA work (Parin et al. 2012).
		core::scoring::ScoreFunctionOP chem_shift_scorefxn = hires_scorefxn_->clone();
		chem_shift_scorefxn->set_weight( rna_chem_shift, CS_weight );
		chem_shift_scorefxn_ = chem_shift_scorefxn;
	}
}

////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::apply( pose::Pose & pose ){

	initialize( pose );

	pose::Pose start_pose = pose;

	monte_carlo_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo( pose, *denovo_scorefxn_, options_->temperature() ) );
	setup_monte_carlo_cycles( pose );

	Size max_tries( 1 );
	if ( options_->filter_lores_base_pairs() || options_->filter_chain_closure() )  max_tries = 10;

	for ( Size ntries = 1; ntries <= max_tries; ++ntries ) {

		time_t pdb_start_time = time(NULL);

		if ( ntries > 1 ) TR << TR.Red << "Did not pass filters. Trying the model again: trial " << ntries << " out of " << max_tries << TR.Reset << std::endl;

		pose = start_pose;

		if ( !refine_pose_ ) rna_structure_parameters_->setup_fold_tree_and_jumps_and_variants( pose );

		rna_structure_parameters_->setup_base_pair_constraints( pose ); // needs to happen after setting cutpoint variants, etc.
		constraint_set_ = pose.constraint_set()->clone();
		rna_chunk_library_->initialize_random_chunks( pose ); //actually not random if only one chunk in each region.

		if ( refine_pose_ ) core::pose::copydofs::copy_dofs_match_atom_names( pose, start_pose );

		if ( options_->close_loops_after_each_move() ) {
			rna_loop_closer_->apply( pose );
			do_close_loops_ = true;
		} else {
			do_close_loops_ = false;
		}

		if ( options_->dump_pdb() ) dump_pdb( pose, "start.pdb" );

		if ( !refine_pose_ ) do_random_moves( pose );

		if ( options_->dump_pdb() ) dump_pdb( pose, "random.pdb" );
		monte_carlo_->reset( pose );

		if ( options_->verbose() ) TR << "Beginning main loop... " << std::endl;

		frag_size_ = 3;

		bool found_solution( true );
		for ( Size r = 1; r <= rounds_; r++ ) {

			if ( options_->verbose() ) TR << TR.Blue << "Beginning round " << r << " of " << rounds_ << TR.Reset << std::endl;

			if ( r == rounds_ && options_->close_loops() ) do_close_loops_ = true;

			//Keep score function coarse for early rounds.
			update_denovo_scorefxn_weights( r );

			monte_carlo_->score_function( *working_denovo_scorefxn_ );

			pose = monte_carlo_->lowest_score_pose();

			// Introduce constraints in stages.
			update_pose_constraints( r, pose );
			monte_carlo_->reset( pose );

			// Finer and finer fragments
			update_frag_size( r );

			// finer rigid body moves
			setup_rigid_body_mover( pose, r ); // needs to happen after fold_tree is decided...

			//////////////////////
			// This is it ... do the loop.
			//////////////////////
			for ( Size i = 1; i <= monte_carlo_cycles_ / rounds_; ++i ) {
				// Make this generic fragment/jump multimover next?
				RNA_move_trial( pose );
			}

			if ( get_native_pose() ) {
				Real const rmsd = core::scoring::all_atom_rmsd( *get_native_pose(), pose );
				if ( options_->verbose() ) TR << "All atom rmsd: " << rmsd << std::endl;
			}

			monte_carlo_->recover_low( pose );
			monte_carlo_->show_counters();
			monte_carlo_->reset_counters();

			if ( r == 2 || rounds_ == 1 ) { //special 'early' stage
				lores_score_early_ = (*working_denovo_scorefxn_)( pose );
				if ( options_->filter_lores_base_pairs_early() ) {
					bool const base_pairs_OK = rna_structure_parameters_->check_base_pairs( pose );
					if ( options_->verbose() ) TR << "Checking base pairs early! Result: " << base_pairs_OK << std::endl;
					if ( !base_pairs_OK ) {
						found_solution = false;
						break;
					}
				}
			}

			if ( r == rounds_ / 2 || rounds_ == 1 ) { // halfway point
				if ( options_->filter_chain_closure_halfway() ) {
					Real const filter_chain_closure_distance_halfway = 2 * options_->filter_chain_closure_distance();
					bool const rna_loops_OK = rna_loop_closer_->check_closure( pose, filter_chain_closure_distance_halfway );
					if ( options_->verbose() ) TR << "Checking loop closure with tolerance of " << filter_chain_closure_distance_halfway << " Angstroms! Result: " << rna_loops_OK << std::endl;
					if ( !rna_loops_OK ) {
						found_solution = false;
						break;
					}
				}
			}
		}

		time_t pdb_end_time = time(NULL);
		if ( options_->verbose() ) TR << "Finished fragment assembly of " << out_file_tag_ << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

		if ( !found_solution ) { // Just try again if early exit from above
			if ( ntries == max_tries ) pose = monte_carlo_->lowest_score_pose();
			continue;
		}
		pose = monte_carlo_->lowest_score_pose();

		if ( options_->close_loops() ) rna_loop_closer_->apply( pose, rna_structure_parameters_->connections() );

		// A bunch of filters
		if ( options_->filter_chain_closure() ) {
			if ( !rna_loop_closer_->check_closure( pose, options_->filter_chain_closure_distance() ) ) {
				if ( options_->verbose() ) TR << "Failed chain closure filter." << std::endl;
				continue;
			}
		}

		if ( options_->filter_lores_base_pairs() ) {
			if ( !rna_structure_parameters_->check_base_pairs( pose ) ) {
				if ( options_->verbose() ) TR << "Failed base pairing filter." << std::endl;
				continue;
			}
		}

		lores_score_final_ = (*working_denovo_scorefxn_)( pose );
		if ( options_->autofilter() ) {
			if ( !check_score_filter( lores_score_final_, all_lores_score_final_ ) ) {
				if ( options_->verbose() ) TR << "Failed score filter." << std::endl;
				continue;
			}
		}
		break; //Pass all the filters, early exit
	} // ++ntries <= max_tries

	// Get the full strength constraint back
	update_denovo_scorefxn_weights( rounds_ );
	update_pose_constraints( rounds_, pose );
	if ( options_->verbose() ) working_denovo_scorefxn_->show( std::cout, pose );
	final_scorefxn_ = working_denovo_scorefxn_;

	lores_pose_ = pose.clone();

	if ( options_->minimize_structure() ) {
		rna_minimizer_->set_allow_insert( rna_structure_parameters_->allow_insert() );
		rna_minimizer_->apply( pose );
		if ( options_->close_loops() ) rna_loop_closer_->apply( pose, rna_structure_parameters_->connections() );
		final_scorefxn_ = hires_scorefxn_;
	}

	if ( options_->use_chem_shift_data() ) apply_chem_shift_data( pose );

	if ( options_->relax_structure() ) rna_relaxer_->apply( pose );

	if ( options_->allow_bulge() ) {
		//Identify and virtual the bulge residues.
		/*Size const num_res_virtualized =*/
		virtualize_bulges(
			pose, options_->allowed_bulge_res(), final_scorefxn_, out_file_tag_,
			true /*allow_pre_virtualize*/, options_->allow_consecutive_bulges(),
			true /*verbose*/
		);
	}

	final_score( pose ); // may include rna_chem_map score here.

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details
// Totally ad hoc recipe for default # cycles, if not user-specified.
void
RNA_FragmentMonteCarlo::setup_monte_carlo_cycles( core::pose::Pose const & pose ){

	if ( options_->user_defined_cycles() || options_->monte_carlo_cycles() > 0 ) {
		monte_carlo_cycles_ = options_->monte_carlo_cycles();
		if ( options_->verbose() ) TR << "Using user-defined " << monte_carlo_cycles_ << " cycles in de novo modeling." << std::endl;
		return;
	}
	// figure out rough number of degrees of freedom.

	// first count up number of residues with allow_insert.
	Size const nres_move = get_moving_res( pose, rna_structure_parameters_->allow_insert() ).size();
	if ( options_->verbose() ) TR << "Number of moving residues: " << nres_move << std::endl;

	Size const nchunks = rna_chunk_library_->num_moving_chunks();
	if ( options_->verbose() ) TR << "Number of moving chunks: " << nchunks << std::endl;

	// then count up rigid bodies that need to be docked.
	Size nbody_move = protocols::farna::get_rigid_body_jumps( pose ).size();
	if ( nbody_move > 1 ) nbody_move--; // first rigid body does not move, by convention.
	if ( nbody_move > 0 ) TR << "Number of moving bodies: " << nbody_move << std::endl;

	monte_carlo_cycles_ = 2000 * ( nres_move + nchunks ) + 20000 * nbody_move;

	if ( monte_carlo_cycles_ > monte_carlo_cycles_max_default_ ) {
		monte_carlo_cycles_ = monte_carlo_cycles_max_default_;
		if ( options_->verbose() ) TR << "Using maximum default Monte Carlo cycles: " <<  monte_carlo_cycles_ << ". Use -cycles option to change this." << std::endl;
	}

	if ( options_->verbose() ) TR << "Using " << monte_carlo_cycles_ << " cycles in de novo modeling." << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::setup_rigid_body_mover( pose::Pose const & pose, Size const r ){

	core::kinematics::MoveMap movemap;
	movemap.set_jump( false );

	bool rigid_body_moves = protocols::farna::let_rigid_body_jumps_move( movemap, pose, options_->move_first_rigid_body() );

	if ( !rigid_body_moves ) return;

	//Keep moves coarse for early rounds. For the last 1/4 of modeling, plateau to the finest moves.
	Real suppress ( 3.0 / 4.0 * r / rounds_ );
	if ( suppress > 1.0 ) suppress = 1.0;

	Real const rot_mag_init( 10.0 ),   rot_mag_final( 0.2 );
	Real const trans_mag_init( 5.0 ), trans_mag_final( 0.1 );
	Real const rot_mag   = rot_mag_init   +  (rot_mag_final - rot_mag_init ) * suppress;
	Real const trans_mag = trans_mag_init +  (trans_mag_final - trans_mag_init ) * suppress;

	rigid_body_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( pose, movemap, rot_mag, trans_mag, protocols::rigid::partner_upstream /*because virtual anchor should be root*/ ) );
	jump_change_frequency_ = 0.5; /* up from default of 0.1*/

	TR << " rot_mag: " << rot_mag << "    trans_mag: " << trans_mag << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::final_score( core::pose::Pose & pose )
{
	core::scoring::ScoreFunctionOP final_scorefxn_plus_slow_terms = final_scorefxn_->clone();
	if ( scoring::rna::nonconst_rna_scoring_info_from_pose( pose ).rna_data_info().rna_reactivities().size() > 0 ) {
		final_scorefxn_plus_slow_terms->set_weight( core::scoring::rna_chem_map, 1.0 );
	}
	( *final_scorefxn_plus_slow_terms )( pose );
	if ( final_scorefxn_plus_slow_terms->has_nonzero_weight( core::scoring::rna_chem_map ) && options_->verbose() ) final_scorefxn_plus_slow_terms->show( pose );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::do_random_moves( core::pose::Pose & pose ) {

	rna_chunk_library_->check_fold_tree_OK( pose );
	rna_chunk_library_->initialize_random_chunks( pose );

	if ( options_->dump_pdb() ) pose.dump_pdb( "add_chunks.pdb" );

	Size const heat_cycles = 3 * pose.total_residue();
	TR << "Heating up... " << std::endl;

	for ( Size i = 1; i <= heat_cycles; i++ ) {
		rna_fragment_mover_->random_fragment_insertion( pose, 1 /*frag_size*/ );
	}

	if ( options_->dump_pdb() )  pose.dump_pdb( "random_moves1.pdb" );

	rna_chunk_library_->initialize_random_chunks( pose, options_->dump_pdb() );

	if ( options_->dump_pdb() )  pose.dump_pdb( "random_moves2.pdb" );

	translate_virtual_anchor_to_first_rigid_body( pose ); //useful for graphics viewing & final output

	if ( options_->dump_pdb() )  pose.dump_pdb( "random_moves3.pdb" );

	randomize_rigid_body_orientations( pose );

	if ( options_->dump_pdb() )  pose.dump_pdb( "random_moves4.pdb" );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::randomize_rigid_body_orientations( pose::Pose & pose ){

	using namespace protocols::rigid;
	using namespace protocols::farna;
	using namespace kinematics;

	utility::vector1< Size > const rigid_body_jumps = get_rigid_body_jumps( pose );
	Size const found_jumps = rigid_body_jumps.size();
	if ( found_jumps <= 1 )  return; // nothing to rotate/translate relative to another object.

	// translation to first, fixed rigid body.
	Vector first_rigid_body_position = pose.jump( rigid_body_jumps[ 1 ] ).get_translation();

	for ( Size n = 2; n <= rigid_body_jumps.size(); n++ ) {
		Size const i = rigid_body_jumps[ n ];

		// randomize orientation.
		RigidBodyRandomizeMover rigid_body_randomize_mover( pose, i, partner_upstream );
		rigid_body_randomize_mover.apply( pose );

		// randomize translation.
		// how far out should we push this segment?
		// For now, hard-wire a value, but later may want to take into account radius of gyration of the chunk.
		Jump jump = pose.jump( i );
		jump.set_translation( first_rigid_body_position );
		pose.set_jump( i, jump ); // move to 'origin' -- position of first rigid body.

		Real const translation_magnitude( 20.0 );
		RigidBodyPerturbMover rigid_body_perturb_mover( i, 0.0 /*rot_mag_in*/, translation_magnitude, partner_upstream );
		rigid_body_perturb_mover.apply( pose );
	}


}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::update_denovo_scorefxn_weights( Size const r )
{
	using namespace core::scoring;

	Real const rna_base_axis_final_weight         = denovo_scorefxn_->get_weight( rna_base_axis );
	Real const rna_base_stagger_final_weight      = denovo_scorefxn_->get_weight( rna_base_stagger );
	Real const rna_base_stack_axis_final_weight   = denovo_scorefxn_->get_weight( rna_base_stack_axis );
	Real const linear_chainbreak_final_weight     = denovo_scorefxn_->get_weight( linear_chainbreak );
	Real const chainbreak_final_weight            = denovo_scorefxn_->get_weight( chainbreak );
	Real const atom_pair_constraint_final_weight  = denovo_scorefxn_->get_weight( atom_pair_constraint );
	Real const coordinate_constraint_final_weight = denovo_scorefxn_->get_weight( coordinate_constraint );
	Real const rna_chem_map_lores_final_weight    = denovo_scorefxn_->get_weight( rna_chem_map_lores );

	//Keep score function coarse for early rounds.
	// Real const suppress  = (r - 1.0) / (rounds - 1.0);
	Real const suppress  = static_cast<Real>( r ) / rounds_;

	working_denovo_scorefxn_->set_weight( rna_base_axis,      suppress*rna_base_axis_final_weight  );
	working_denovo_scorefxn_->set_weight( rna_base_stagger,   suppress*rna_base_stagger_final_weight  );
	if ( options_->titrate_stack_bonus() ) working_denovo_scorefxn_->set_weight( rna_base_stack_axis,suppress*rna_base_stack_axis_final_weight  );
	working_denovo_scorefxn_->set_weight( atom_pair_constraint,  suppress*atom_pair_constraint_final_weight  );
	working_denovo_scorefxn_->set_weight( coordinate_constraint,  suppress*coordinate_constraint_final_weight  );
	working_denovo_scorefxn_->set_weight( rna_chem_map_lores,   suppress*rna_chem_map_lores_final_weight  );

	// keep chainbreak extra low for early rounds... seems to be important for rigid body modeler.
	//Real suppress_chainbreak  = ( r - ( rounds_/3.0 ) )/ ( static_cast<Real>( rounds_ ) - ( rounds_ / 3.0 ) );
	//Real const suppress_chainbreak_min = 1 / static_cast< Real >( rounds_ );
	//if ( suppress_chainbreak < suppress_chainbreak_min ) suppress_chainbreak = suppress_chainbreak_min;

	working_denovo_scorefxn_->set_weight( linear_chainbreak,  suppress * linear_chainbreak_final_weight  );
	working_denovo_scorefxn_->set_weight( chainbreak,  suppress * chainbreak_final_weight  );

}


////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_FragmentMonteCarlo::figure_out_constraint_separation_cutoff( Size const r, Size const max_dist )
{

	//Keep score function coarse for early rounds.
	Real const suppress  = 5.0 / 3.0 * r / rounds_;

	Size const separation_cutoff ( static_cast<Size>( suppress * max_dist ) + 2 );
	if ( separation_cutoff > max_dist ) return max_dist;

	return separation_cutoff;
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::update_pose_constraints( Size const r, core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;

	if ( !options_->staged_constraints() ) return;

	if ( !constraint_set_ ) return;

	static core::kinematics::ShortestPathInFoldTree shortest_path_in_fold_tree( pose.fold_tree() );
	Size const separation_cutoff = figure_out_constraint_separation_cutoff( r, shortest_path_in_fold_tree.max_dist() );
	TR << "ROUND " << r << " out of " << rounds_ << std::endl;
	TR << "FOLD_TREE CURRENT SEPARATION CUTOFF " << separation_cutoff << " out of " << shortest_path_in_fold_tree.max_dist() << std::endl;
	// Fang: apparently separation_cutoff can be smaller than shortest_path_in_fold_tree.dist( i , j )
	// Not sure why but do an early exit hack here
	if ( separation_cutoff >= shortest_path_in_fold_tree.max_dist() ) {
		pose.constraint_set( constraint_set_ );
	} else {
		ConstraintCOPs csts( constraint_set_->get_all_constraints() );
		ConstraintSetOP cst_set_new( new scoring::constraints::ConstraintSet );
		for ( Size n = 1; n <= csts.size(); n++ ) {
			ConstraintCOP const & cst( csts[n] );
			if ( cst->natoms() == 2 )  { // currently only defined for pairwise distance constraints.
				Size const i = cst->atom( 1 ).rsd();
				Size const j = cst->atom( 2 ).rsd();
				Size const dist( shortest_path_in_fold_tree.dist( i , j ) );
				if ( dist  > separation_cutoff ) continue;
			}
			cst_set_new->add_constraint( cst );
		}
		pose.constraint_set( cst_set_new );
	}

	TR << "NUM CONSTRAINTS " << pose.constraint_set()->get_all_constraints().size() << " out of " <<
		constraint_set_->get_all_constraints().size() << std::endl;

}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::update_frag_size( Size const r )
{
	Real const frag_size_old = frag_size_;
	frag_size_ = 3;
	if ( r > 1.0 * ( rounds_ / 3.0 ) ) frag_size_ = 2;
	if ( r > 2.0 * ( rounds_ / 3.0 ) ) frag_size_ = 1;
	if ( frag_size_ != frag_size_old ) TR << "Fragment size: " << frag_size_ << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////
/// @details
/// There are now two kinds of insertions --
/// (1) fragment insertions for, e.g., contiguous 3-mers
///   and
/// (2) "chunk insertions", which change out whole loops, motifs, or
///     junctions based on previous models stored in silent files
void
RNA_FragmentMonteCarlo::RNA_move_trial( pose::Pose & pose ) {



	if  ( numeric::random::rg().uniform() < jump_change_frequency_ )  {
		//Following returns early if there are no jumps.
		random_jump_trial( pose );
	} else {

		bool did_a_trial( false );
		if ( numeric::random::rg().uniform() < chunk_coverage_ ) {
			did_a_trial = random_chunk_trial( pose );
		}

		if ( !did_a_trial ) {
			random_fragment_trial( pose );
		}
	}


}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::random_jump_trial( pose::Pose & pose ) {

	bool success( false );
	std::string move_type( "" );

	if ( rigid_body_mover_ &&  numeric::random::rg().uniform() < 0.8 /*totally arbitrary*/ ) {
		rigid_body_mover_->apply( pose );
		success = true; /* rigid body mover is from docking  */
		move_type = "rigid_body";
	} else {
		success = rna_structure_parameters_->random_jump_change( pose );
		move_type = "jump_change";
	}

	if ( !success ) return;

	if ( do_close_loops_ )  rna_loop_closer_->apply( pose );

	monte_carlo_->boltzmann( pose, move_type );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::random_fragment_trial( pose::Pose & pose ) {

	rna_fragment_mover_->random_fragment_insertion( pose, frag_size_ );
	if ( do_close_loops_ ) rna_loop_closer_->apply( pose );

	monte_carlo_->boltzmann( pose, "frag" + SS(frag_size_) );

}

////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_FragmentMonteCarlo::random_chunk_trial( pose::Pose & pose ) {

	bool const did_an_insertion = rna_chunk_library_->random_chunk_insertion( pose );
	if ( !did_an_insertion ) return false;

	if ( do_close_loops_ ) rna_loop_closer_->apply( pose );

	monte_carlo_->boltzmann( pose, "chunk" );

	return true /*did an insertion*/;

}

///////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_FragmentMonteCarlo::check_score_filter( Real const lores_score, std::list< Real > & all_lores_score ){

	all_lores_score.push_back( lores_score );

	all_lores_score.sort(); // nice -- can do this with a list!

	// note that if options_->autofilter_score_quantile_ = 0.20, the first decoy will be 'passed' for free.
	Size const n = all_lores_score.size();
	Size const cutoff_index = static_cast< Size >( n * options_->autofilter_score_quantile() ) + 1;

	// the one pain with lists -- need to iterate through to find the element corresponding to the quantile score.
	Real all_lores_score_cutoff = all_lores_score.front();
	Size i( 1 );
	for ( std::list< Real >::const_iterator iter = all_lores_score.begin(), end = all_lores_score.end(); iter != end; ++iter, i++ ) {
		if ( i == cutoff_index ) all_lores_score_cutoff = *iter;
	}

	TR << "Comparing current lores score " << lores_score << " to automatically determined cutoff: " << all_lores_score_cutoff << " based on " << options_->autofilter_score_quantile() << " quantile from "  << n << " models so far" << std::endl;
	return ( lores_score <= all_lores_score_cutoff );
}

///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarlo::apply_chem_shift_data( core::pose::Pose & pose ){

	using namespace core::scoring;
	using namespace core::io::pdb;

	runtime_assert( options_->use_chem_shift_data() );

	if ( options_->minimize_structure() ) {
		rna_minimizer_->set_score_function( chem_shift_scorefxn_ ); //use the chem_shift_scorefxn_
		rna_minimizer_->apply( pose );
		rna_minimizer_->set_score_function( hires_scorefxn_ ); //set back the original scorefxn.
		if ( options_->close_loops() ) rna_loop_closer_->apply( pose, rna_structure_parameters_->connections() );
		final_scorefxn_ = chem_shift_scorefxn_;
	}

	(*chem_shift_scorefxn_)( pose );
}

void
RNA_FragmentMonteCarlo::show(std::ostream & output) const
{
	Mover::show(output);
	output <<  "\nRounds:                        " << rounds_ <<
		"\nMonte Carlo cycles:            " << monte_carlo_cycles_ <<
		"\nMC cycle max default:          " << monte_carlo_cycles_max_default_ <<
		"\nFragment size:                 " << frag_size_ <<
		"\nDo close loops:                " << (do_close_loops_ ? "True" : "False") <<
		"\nJump change frequency:         " << jump_change_frequency_ <<
		"\nLores score early:             " << (lores_score_early_) <<
		"\nLores score final:             " << (lores_score_final_) <<
		"\nChunk coverage:                " << chunk_coverage_;
}

} //farna
} //protocols
