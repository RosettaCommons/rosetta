// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/movers/RNA_DeNovoMasterMover.cc
/// @brief  Fragment insertions, jump moves, docking, 'chunk' moves, and RNA-protein docking
/// @author Rhiju Das, rhiju@stanford.edu
/// @author Kalli Kappel, kappel@stanford.edu

#include <protocols/rna/denovo/movers/RNA_DeNovoMasterMover.hh>
#include <core/import_pose/options/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/rna/denovo/movers/RNA_DeNovoMasterMover.hh>
#include <protocols/rna/denovo/movers/RNA_FragmentMover.hh>
#include <core/fragment/rna/RNA_Fragments.hh>
#include <core/fragment/rna/FullAtomRNA_Fragments.hh>
#include <core/import_pose/RNA_BasePairHandler.hh>
#include <core/import_pose/libraries/RNA_JumpLibrary.hh>
#include <core/import_pose/libraries/RNA_ChunkLibrary.hh>
#include <core/import_pose/libraries/RNA_LibraryManager.hh>
#include <core/import_pose/RNA_JumpMover.hh>
#include <protocols/rna/denovo/movers/RNA_HelixMover.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rna/movers/RNA_LoopCloser.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/format.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility>

static basic::Tracer TR( "protocols.rna.denovo.movers.RNA_DeNovoMasterMover" );

using namespace core;
using namespace core::pose::rna;
using namespace ObjexxFCL::format; // AUTO USING NS
using namespace core::fragment::rna;
using namespace core::import_pose;
using namespace core::import_pose::libraries;
using utility::vector1;

/////////////////////////////////////////////////////////////////////////
/// @detailed
///   RNA_DeNovoMasterMover (factored out of RNA_FragmentMonteCarlo)
///
///    Decides whether or not to run fragment insertions, RNA jump moves
///           docking, 'chunk' moves, etc.
///
///    This class should *not* have to know details of MonteCarlo or
///         simulated tempering scheme -- have that handled outside.
///
///  TODO:
///   Better unify RNA & RNP docking movers. [perhaps set those up at beginning, and
///          only update rot/trans mags in here?]
///
/////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

//Constructor
RNA_DeNovoMasterMover::RNA_DeNovoMasterMover( core::import_pose::options::RNA_FragmentMonteCarloOptionsCOP options,
	core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::import_pose::RNA_BasePairHandlerCOP rna_base_pair_handler,
	protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer,
	libraries::RNA_ChunkLibraryOP rna_chunk_library ):
	options_(std::move( options )),
	frag_size_( 3 ),
	jump_change_frequency_( 0.1 ), //  maybe updated based on options, or if rigid-body sampling
	dock_into_density_freq_( 0.2 ),
	close_loops_( false ),
	do_rnp_docking_( false ),
	move_type_( "" ),
	rna_chunk_library_(std::move( rna_chunk_library )),
	rna_loop_closer_(std::move( rna_loop_closer ))
{
	RNA_Fragments const & all_rna_fragments_( RNA_LibraryManager::get_instance()->rna_fragment_library( options_->all_rna_fragments_file() ) );
	rna_fragment_mover_ = RNA_FragmentMoverOP( new RNA_FragmentMover( all_rna_fragments_, atom_level_domain_map, options_->symm_hack_arity() ) );

	rna_jump_mover_ = RNA_JumpMoverOP( new RNA_JumpMover( RNA_LibraryManager::get_instance()->rna_jump_library_cop( options_->jump_library_file() ),
		atom_level_domain_map ) );
	rna_jump_mover_->set_rna_pairing_list ( rna_base_pair_handler->rna_pairing_list() );
	rna_jump_mover_->set_chain_connections( rna_base_pair_handler->chain_connections() );
	jump_change_frequency_ = options_->jump_change_frequency(); // might get overwritten if rigid_body movement, tested later.
}

//Destructor
RNA_DeNovoMasterMover::~RNA_DeNovoMasterMover() = default;

void
RNA_DeNovoMasterMover::apply( core::pose::Pose & pose ) {
	apply( pose, 1 );
}

void
RNA_DeNovoMasterMover::apply( core::pose::Pose & pose,
	Size const & cycle_number )
{
	do_move_trial( cycle_number, pose );
}

////////////////////////////////////////////////////////////////////////////////////////
/// @details
///    Do a move.
///    RNA moves are fragment/jump/chunks.
///    Every 10th iteration, try a RNA/protein docking move (K. Kappel)
///
void
RNA_DeNovoMasterMover::do_move_trial( Size const & i, pose::Pose & pose )
{

	move_type_ = ""; // if left blank, no move was done.

	// If RNA/protein, do docking every 10th move) - remove once we no longer need
	// rna_protein_docking_legacy
	if ( do_rnp_docking_ && (i % 10 /* old default for rna_protein_docking_freq*/ == 0 )
			&& options_->rna_protein_docking_legacy() ) {

		rnp_docking_trial( pose );

	} else {
		// regular RNA move (fragment, chunk, or jump)
		RNA_move_trial( pose );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
/// @details
/// There are three kinds of insertions --
/// (1) fragment insertions for, e.g., contiguous 3-mers
/// (2) jump changes (tweaking base pairs)
///   and
/// (3) "chunk insertions", which change out whole loops, motifs, or
///     junctions based on previous models stored in silent files
void
RNA_DeNovoMasterMover::RNA_move_trial( pose::Pose & pose ) {


	/// if density map was provided:
	/// ~20% (totally arbitrary) of moves should be rigid body docking into density
	/// (LEGACY)
	if ( options_->dock_into_density() && options_->dock_into_density_legacy() &&
			numeric::random::rg().uniform() < dock_into_density_freq_ ) {

		dock_into_density_trial( pose );

	} else if ( numeric::random::rg().uniform() < jump_change_frequency_ )  {
		//Following returns early if there are no jumps.
		//does rigid body docking as well as random jump insertions
		random_jump_trial( pose );

	} else {

		bool did_a_trial( false );
		if ( numeric::random::rg().uniform() < rna_chunk_library_->chunk_coverage() ) {
			did_a_trial = random_chunk_trial( pose );
		}

		if ( !did_a_trial ) {
			random_fragment_trial( pose );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::random_jump_trial( pose::Pose & pose ) {

	bool success( false );
	std::string move_type( "" );

	if ( rna_helix_mover_ && numeric::random::rg().uniform() < 0.25 /*made up*/ ) {
		rna_helix_mover_->apply( pose );
		success = true; // wait actually we should check in helix mover whether stuff actually moved
		move_type = "helix_move";
	} else if ( rigid_body_mover_ &&  numeric::random::rg().uniform() < 0.8 /*totally arbitrary*/ && options_->docking() ) {
		rigid_body_mover_->apply( pose );
		success = true; /* rigid body mover is from docking  */
		move_type = "rigid_body";
	} else {
		success = rna_jump_mover_->random_jump_change( pose );
		move_type = "jump_change";
	}

	if ( !success ) return;

	if ( close_loops_ )  rna_loop_closer_->apply( pose );

	move_type_ = move_type;
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::dock_into_density_trial( pose::Pose & pose ) {

	std::string move_type( "dock_density" );
	dock_into_density_mover_->apply( pose );
	move_type_ = move_type;

}
////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::rnp_docking_trial( pose::Pose & pose ) {

	// clear energies after applying a new fold tree to the pose
	// to prevent any scoring craziness

	if ( !rnp_docking_mover_ ) return;
	std::string move_type( "rnp_dock" );

	// initial fold tree
	core::kinematics::FoldTree const ft_init = pose.fold_tree();
	// apply the docking fold tree
	// jumps only between RNA and protein
	pose.fold_tree( rnp_docking_ft_ );
	// don't need to clear here, because we won't score until after fold tree
	// is returned and we'll clear energies then

	rnp_docking_mover_->apply( pose );

	// Need to return the fold tree before scoring
	// so that linear_chainbreak gets computed correctly...
	// there must be a better way to do this, but we need
	// the fold tree set up so protein and RNA chains are
	// separated by jumps (with no jumps within RNA or
	// protein chains) so that the chains get docked as a
	// whole
	// could probably set up the docking mover in a different
	// way... but this works for now

	pose.fold_tree( ft_init );
	pose.energies().clear(); // clear the energies so that pose is scored correctly next time

	move_type_ = move_type;

}
////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::random_fragment_trial( pose::Pose & pose ) {

	rna_fragment_mover_->random_fragment_insertion( pose, frag_size_ );
	if ( close_loops_ ) rna_loop_closer_->apply( pose );

	move_type_ = "frag" + SS(frag_size_);

}

////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_DeNovoMasterMover::random_chunk_trial( pose::Pose & pose ) {

	bool const did_an_insertion = rna_chunk_library_->random_chunk_insertion( pose );
	if ( !did_an_insertion ) return false;

	if ( close_loops_ ) rna_loop_closer_->apply( pose );

	move_type_ = "chunk";

	return true /*did an insertion*/;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::setup_rigid_body_mover( pose::Pose const & pose,
	Real const & rot_mag,
	Real const & trans_mag )
{

	core::kinematics::MoveMap movemap;
	movemap.set_jump( false );

	bool move_first_rigid_body = false;
	if ( options_->move_first_rigid_body() || options_->dock_into_density() ) {
		move_first_rigid_body = true;
	}

	bool allow_single_rigid_body = false;
	if ( options_->dock_into_density() ) allow_single_rigid_body = true;

	// this updates the movemap
	bool rigid_body_moves = let_rigid_body_jumps_move( movemap, pose, move_first_rigid_body, allow_single_rigid_body /* for moving absolute coordinates*/ );

	if ( !rigid_body_moves ) return;

	rigid_body_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( pose, movemap, rot_mag, trans_mag, protocols::rigid::partner_upstream /*because virtual anchor should be root*/ ) );
	if ( options_->rna_protein_docking_freq() !=  0.4 /* default */ ) {
		// if the user specified rna protein docking frequency
		// then use this to set the jump change frequency, i.e. how
		// often we'll try to do rigid body moves
		jump_change_frequency_ = options_->rna_protein_docking_freq() / 0.8;
	} else {
		jump_change_frequency_ = 0.5; /* up from default of 0.1*/
	}

	TR << " rot_mag: " << rot_mag << "    trans_mag: " << trans_mag << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::set_helix_mover_magnitude(
	core::Real const & rot_mag,
	core::Real const & trans_mag )
{
	if ( !rna_helix_mover_ ) return;

	rna_helix_mover_->set_rot_magnitude( rot_mag );
	rna_helix_mover_->set_trans_magnitude( trans_mag );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Only used for legacy RNP setup, will remove once all testing is complete
void
RNA_DeNovoMasterMover::setup_dock_into_density_mover( pose::Pose const & pose,
	Real const & rot_mag,
	Real const & trans_mag )
{
	if ( !options_->dock_into_density() ) return;

	core::kinematics::MoveMap movemap;
	movemap.set_jump( false );

	//Figure out the jump that we want to sample
	// Fold tree better already be set up properly!!!
	for ( core::Size i=1; i<=pose.fold_tree().num_jump(); ++i ) {
		Size upstream_res = pose.fold_tree().upstream_jump_residue( i );
		if ( pose.residue( upstream_res ).aa() == core::chemical::aa_vrt ) {
			// this is the jump from the virtual residue
			movemap.set_jump( i, true );
			break;
		}
	}

	dock_into_density_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( pose, movemap, rot_mag, trans_mag, protocols::rigid::partner_upstream /*because virtual anchor should be root*/ ) );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Kalli -- this copies code from setup_rigid_body_mover -- can we unify? -- rhiju, 2017
/// yes, unified! but still need to keep this here for legacy RNP setup. will remove later
/// once all testing is complete  - Kalli, late 2017
void
RNA_DeNovoMasterMover::setup_rna_protein_docking_mover( pose::Pose const & pose,
	Real const & rot_mag,
	Real const & trans_mag )
{
	do_rnp_docking_ = false;

	core::kinematics::MoveMap movemap;
	movemap.set_jump( false );

	if ( rnp_docking_ft_.num_jump() < 1 ) return;

	do_rnp_docking_ = true;

	for ( core::Size i =1; i<=rnp_docking_ft_.num_jump(); ++i ) {
		// Double check that there's RNA on one side and protein on the other
		Size up_res = rnp_docking_ft_.upstream_jump_residue( i );
		Size down_res = rnp_docking_ft_.downstream_jump_residue( i );
		if ( (pose.residue( up_res ).is_RNA() && pose.residue( down_res ).is_protein() ) ||
				(pose.residue( up_res ).is_protein() && pose.residue( down_res ).is_RNA()) ) {
			movemap.set_jump( i, true );
			TR << "Set up RNP docking jump between residue " << up_res << " and " << down_res << std::endl;
		}
	}
	rnp_docking_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( pose, movemap, rot_mag, trans_mag, protocols::rigid::partner_upstream ) );
}

////////////////////////////////////////////////////////////////////////////////////////
/// @details
/// Initial 'heating' of model when starting from scratch -- a bunch of
///  fragment insertions, jumps, chunks, etc. with no energy function checks.
void
RNA_DeNovoMasterMover::do_random_moves( core::pose::Pose & pose,
	Size const & monte_carlo_cycles,
	bool const & check_num_rna_res /* = false*/ ) {

	rna_chunk_library_->check_fold_tree_OK( pose );
	rna_chunk_library_->initialize_random_chunks( pose );

	if ( options_->dump_pdb() ) pose.dump_pdb( "add_chunks.pdb" );

	translate_virtual_anchor_to_first_rigid_body( pose ); //useful for graphics viewing & final output

	Size nres( 0 );
	if ( check_num_rna_res ) {
		Size num_RNA_res = 0;
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.residue( i ).is_RNA() ) num_RNA_res += 1;
		}
		nres = num_RNA_res;
	} else {
		nres = pose.size();
	}


	Size const heat_cycles = std::min( 3 * nres, monte_carlo_cycles );
	TR << "Heating up... " << heat_cycles << " cycles." << std::endl;

	for ( Size i = 1; i <= heat_cycles; i++ ) {
		rna_fragment_mover_->random_fragment_insertion( pose, 1 /*frag_size*/ );
	}

	rna_chunk_library_->initialize_random_chunks( pose );

	translate_virtual_anchor_to_first_rigid_body( pose ); //useful for graphics viewing & final output

	randomize_rigid_body_orientations( pose );

	if ( options_->dump_pdb() ) pose.dump_pdb( "random_moves.pdb" );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::randomize_rigid_body_orientations( pose::Pose & pose ){

	using namespace protocols::rigid;
	using namespace protocols::rna::denovo;
	using namespace kinematics;

	vector1< Size > const rigid_body_jumps = get_rigid_body_jumps( pose );
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
// make this a separate mover?
void
RNA_DeNovoMasterMover::search_rigid_body_orientation( pose::Pose & pose ){

	using namespace core::scoring;
	using namespace protocols::rigid;

	// find a good starting configuration - designed initially with RNP density docking in mind
	// for now only for RNPs, but once I fix fold tree issues this could easily
	// generalize to any type of partners

	// get score function
	// by default use a score function that just has VDW and density score terms
	ScoreFunctionOP scorefxn;

	// make option to set this as something else
	scorefxn = ScoreFunctionFactory::create_score_function( "rna/denovo/rna_and_rnp_vdw.wts" );
	if ( options_->model_with_density() ) {
		// turn on elec_dens_fast
		if ( !scorefxn->has_nonzero_weight( elec_dens_fast ) ) scorefxn->set_weight( elec_dens_fast, 10.0 );
	}

	// make monte carlo object

	core::kinematics::FoldTree const ft_init = pose.fold_tree();
	pose.fold_tree( rnp_docking_ft_ );
	protocols::moves::MonteCarloOP monte_carlo( new protocols::moves::MonteCarlo( pose, *scorefxn, 2.0 ) );

	// make rigid body mover
	if ( rnp_docking_ft_.num_jump() < 1 ) return;

	utility::vector1< Size > docking_jumps;
	for ( Size i = 1; i <= rnp_docking_ft_.num_jump(); ++i ) {
		Size up_res = rnp_docking_ft_.upstream_jump_residue( i );
		Size down_res = rnp_docking_ft_.downstream_jump_residue( i );
		if ( ! ((pose.residue( up_res ).is_RNA() && pose.residue( down_res ).is_protein() ) ||
				(pose.residue( up_res ).is_protein() && pose.residue( down_res ).is_RNA())) ) {
			continue;
		}
		docking_jumps.push_back( i );
	}

	core::Real rot_mag( 10.0 );
	core::Real trans_mag( 5.0 );

	// do a bunch of random rigid body moves
	// this only works for systems with one RNP jump, it's set up only for the legacy rnp fold tree stuff
	for ( Size jump_index = 1; jump_index <= docking_jumps.size(); ++jump_index ) {
		Size jump = docking_jumps[ jump_index ];

		// randomize the orientation
		RigidBodyRandomizeMover rigid_body_randomize_mover( pose, jump, partner_upstream );

		// and the do a random perturbation
		core::kinematics::MoveMap movemap;
		movemap.set_jump( jump, true );
		RigidBodyPerturbMoverOP random_perturb_mover = RigidBodyPerturbMoverOP( new RigidBodyPerturbMover( pose, movemap, rot_mag, trans_mag, protocols::rigid::partner_downstream ) );

		core::pose::Pose start_pose = pose;
		for ( Size i = 1; i <= 200; ++i ) {
			pose = start_pose;
			pose.fold_tree( rnp_docking_ft_ );

			rigid_body_randomize_mover.apply( pose );

			random_perturb_mover->apply( pose );
			// need init fold tree for density scoring
			pose.fold_tree( ft_init );
			pose.energies().clear(); // clear the energies so that pose is scored correctly next time

			monte_carlo->boltzmann( pose, "search_rigid" );
		}
	}

	// pick the pose with the best score
	// (we're basically trying to find places with open density)
	monte_carlo->show_counters();
	pose = monte_carlo->lowest_score_pose();
	// make sure it has the correct fold tree
	pose.fold_tree( ft_init );
	pose.energies().clear();


}
////////////////////////////////////////////////////////////////////////////////////////
/// Kalli -- this copies code from randomize_rigid_body_orientations -- can we unify? -- rhiju, 2017
/// yes, unified! but still need to keep this here for legacy RNP setup. will remove later
/// once all testing is complete  - Kalli, late 2017
void
RNA_DeNovoMasterMover::randomize_rnp_rigid_body_orientations( pose::Pose & pose ){

	using namespace protocols::rigid;
	using namespace protocols::rna::denovo;
	using namespace kinematics;

	if ( rnp_docking_ft_.num_jump() < 1 ) return;

	for ( Size i = 1; i <= rnp_docking_ft_.num_jump(); ++i ) {
		// check that upstream and downstream are rna and protein
		// if we are doing density modeling, this shouldn't be true for the last jump
		Size up_res = rnp_docking_ft_.upstream_jump_residue( i );
		Size down_res = rnp_docking_ft_.downstream_jump_residue( i );
		if ( ! ((pose.residue( up_res ).is_RNA() && pose.residue( down_res ).is_protein() ) ||
				(pose.residue( up_res ).is_protein() && pose.residue( down_res ).is_RNA())) ) {
			continue;
		}

		pose::Pose pose_copy = pose;

		// randomize orientation
		RigidBodyRandomizeMover rigid_body_randomize_mover( pose, i );
		rigid_body_randomize_mover.apply( pose );
		RigidBodyRandomizeMover rigid_body_randomize_mover_2( pose, i, partner_upstream );
		rigid_body_randomize_mover_2.apply( pose );
	}
}

////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::setup_rnp_fold_tree( pose::Pose & pose, bool const & refine_pose, bool const & randomize )
{
	// set up the fold tree for docking
	rnp_docking_ft_ = get_rnp_docking_fold_tree( pose, options_->model_with_density() );

	if ( !refine_pose && randomize ) {
		kinematics::FoldTree ft_init = pose.fold_tree(); // save original fold tree
		// apply docking fold tree
		pose.fold_tree( rnp_docking_ft_ );
		pose.energies().clear(); // clear the energies so that pose is scored correctly next time
		// randomize the RNA/protein rigid body orientations
		// IF we're doing RNA/protein docking
		randomize_rnp_rigid_body_orientations( pose );
		// get back original fold tree
		pose.fold_tree( ft_init );
		pose.energies().clear();
		if ( options_->dump_pdb() )  pose.dump_pdb( "random_rnp_orientation.pdb" );
	}
}

} //movers
} //denovo
} //rna
} //protocols
