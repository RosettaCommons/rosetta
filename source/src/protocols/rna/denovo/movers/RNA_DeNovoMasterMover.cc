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
/// @author Kalli Kappel, kkappel@stanford.edu

#include <protocols/rna/denovo/movers/RNA_DeNovoMasterMover.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/rna/denovo/movers/RNA_DeNovoMasterMover.hh>
#include <protocols/rna/denovo/movers/RNA_FragmentMover.hh>
#include <protocols/rna/denovo/fragments/RNA_Fragments.hh>
#include <protocols/rna/denovo/fragments/FullAtomRNA_Fragments.hh>
#include <protocols/rna/denovo/base_pairs/RNA_BasePairHandler.hh>
#include <protocols/rna/denovo/libraries/RNA_JumpLibrary.hh>
#include <protocols/rna/denovo/libraries/RNA_ChunkLibrary.hh>
#include <protocols/rna/denovo/libraries/RNA_LibraryManager.hh>
#include <protocols/rna/denovo/movers/RNA_JumpMover.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rna/movers/RNA_LoopCloser.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/format.hh>
static basic::Tracer TR( "protocols.rna.denovo.movers.RNA_DeNovoMasterMover" );

using namespace core;
using namespace ObjexxFCL::format; // AUTO USING NS
using namespace protocols::rna::denovo::fragments;
using namespace protocols::rna::denovo::libraries;
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
	RNA_DeNovoMasterMover::RNA_DeNovoMasterMover( options::RNA_FragmentMonteCarloOptionsCOP options,
																								protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
																								base_pairs::RNA_BasePairHandlerCOP rna_base_pair_handler,
																								protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer,
																								libraries::RNA_ChunkLibraryOP rna_chunk_library ):
		options_( options ),
		frag_size_( 3 ),
		jump_change_frequency_( 0.1 ), //  maybe updated based on options, or if rigid-body sampling
		close_loops_( false ),
		do_rnp_docking_( false ),
		move_type_( "" ),
		rna_chunk_library_( rna_chunk_library ),
		rna_loop_closer_( rna_loop_closer )
	{
		RNA_Fragments const & all_rna_fragments_( RNA_LibraryManager::get_instance()->rna_fragment_library( options_->all_rna_fragments_file() ) );
		rna_fragment_mover_ = RNA_FragmentMoverOP( new RNA_FragmentMover( all_rna_fragments_, atom_level_domain_map ) );

		rna_jump_mover_ = RNA_JumpMoverOP( new RNA_JumpMover( RNA_LibraryManager::get_instance()->rna_jump_library_cop( options_->jump_library_file() ),
																													atom_level_domain_map ) );
		rna_jump_mover_->set_rna_pairing_list ( rna_base_pair_handler->rna_pairing_list() );
		rna_jump_mover_->set_chain_connections( rna_base_pair_handler->chain_connections() );
		jump_change_frequency_ = options_->jump_change_frequency(); // might get overwritten if rigid_body movement, tested later.
	}

	//Destructor
	RNA_DeNovoMasterMover::~RNA_DeNovoMasterMover()
	{}

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
	// If RNA/protein, do docking every 10th move)
	if ( do_rnp_docking_ && (i % options_->rna_protein_docking_freq() == 0 ) ) {
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

	if  ( numeric::random::rg().uniform() < jump_change_frequency_ )  {
		//Following returns early if there are no jumps.
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

	if ( rigid_body_mover_ &&  numeric::random::rg().uniform() < 0.8 /*totally arbitrary*/ ) {
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
	//ft.show( std::cout );

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

	bool rigid_body_moves = let_rigid_body_jumps_move( movemap, pose, options_->move_first_rigid_body() );

	if ( !rigid_body_moves ) return;

	rigid_body_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( pose, movemap, rot_mag, trans_mag, protocols::rigid::partner_upstream /*because virtual anchor should be root*/ ) );
	jump_change_frequency_ = 0.5; /* up from default of 0.1*/

	TR << " rot_mag: " << rot_mag << "    trans_mag: " << trans_mag << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Kalli -- this copies code from setup_rigid_body_mover -- can we unify? -- rhiju, 2017
void
RNA_DeNovoMasterMover::setup_rna_protein_docking_mover( pose::Pose const & pose,
																							 Real const & rot_mag,
																							 Real const & trans_mag )
{
	do_rnp_docking_ = false;

	core::kinematics::MoveMap movemap;
	movemap.set_jump( false );

	// // figure out the jumps between RNA and protein
	// This doesn't work because there are often jumps within RNA chains, then the whole RNA chain doesn't
	// dock together as a rigid body...
	//  But ultimately this may be a better way to go
	// vector1< Size > rna_protein_jumps;
	// for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {
	//  TR.Debug << "checking jump: " <<  pose.fold_tree().upstream_jump_residue( n ) << " to " <<  pose.fold_tree().downstream_jump_residue( n ) << std::endl;
	//  // if the upstream/downstream jump residues are RNA and protein, then we can move this jump
	//  if (( pose.residue( pose.fold_tree().upstream_jump_residue( n ) ).is_RNA() &&
	//   pose.residue( pose.fold_tree().downstream_jump_residue( n ) ).is_protein() ) ||
	//   (pose.residue( pose.fold_tree().upstream_jump_residue( n ) ).is_protein() &&
	//   pose.residue( pose.fold_tree().upstream_jump_residue( n ) ).is_RNA() )) {
	//
	//    rna_protein_jumps.push_back( n );
	//    std::cout << "Found RNA/protein jump between " << pose.fold_tree().downstream_jump_residue( n ) <<
	//     " and " << pose.fold_tree().upstream_jump_residue( n ) << std::endl;
	//
	//  }
	// }

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
																				Size const & monte_carlo_cycles ) {

	rna_chunk_library_->check_fold_tree_OK( pose );
	rna_chunk_library_->initialize_random_chunks( pose );

	if ( options_->dump_pdb() ) pose.dump_pdb( "add_chunks.pdb" );

	translate_virtual_anchor_to_first_rigid_body( pose ); //useful for graphics viewing & final output

	Size const heat_cycles = std::min( 3 * pose.size(), monte_carlo_cycles );
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
/// Kalli -- this copies code from randomize_rigid_body_orientations -- can we unify? -- rhiju, 2017
void
RNA_DeNovoMasterMover::randomize_rnp_rigid_body_orientations( pose::Pose & pose ){

	using namespace protocols::rigid;
	using namespace protocols::rna::denovo;
	using namespace kinematics;

	if ( rnp_docking_ft_.num_jump() < 1 ) return;

	for ( Size i = 1; i <= rnp_docking_ft_.num_jump(); ++i ) {

		Vector rigid_body_position = pose.jump( i ).get_translation();

		// randomize orientation
		RigidBodyRandomizeMover rigid_body_randomize_mover( pose, i );
		rigid_body_randomize_mover.apply( pose );

		//randomize translation
		Jump jump = pose.jump( i );
		jump.set_translation( rigid_body_position );
		pose.set_jump( i, jump );

		Real const translation_magnitude( 20.0 );
		RigidBodyPerturbMover rigid_body_perturb_mover( i /*jump*/, 0.0 /*rotation*/, translation_magnitude );
		rigid_body_perturb_mover.apply( pose );
	}
}

////////////////////////////////////////////////////////////////////
void
RNA_DeNovoMasterMover::setup_rnp_fold_tree( pose::Pose & pose, bool const & refine_pose )
{
	// set up the fold tree for docking
	rnp_docking_ft_ = get_rnp_docking_fold_tree( pose );
	if ( !refine_pose ) {
		kinematics::FoldTree ft_init = pose.fold_tree(); // save original fold tree
		// apply docking fold tree
		pose.fold_tree( rnp_docking_ft_ );
		pose.energies().clear(); // clear the energies so that pose is scored correctly next time
		// randomize the RNA/protein rigid body orientations
		randomize_rnp_rigid_body_orientations( pose );
		// get back original fold tree
		pose.fold_tree( ft_init );
		pose.energies().clear();
		if ( options_->dump_pdb() )  pose.dump_pdb( "random_rnp_orientation.pdb" );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
kinematics::FoldTree
RNA_DeNovoMasterMover::get_rnp_docking_fold_tree( pose::Pose const & pose ) {

	// This is super simple right now, later update this so that it can take
	// user input docking jumps, that's sort of the next simplest thing that
	// we can do, but eventually can try to be smart looking at chains and
	// distances (a little difficult if there are e.g. 2 RNA chains in a helix
	// that together bind to the protein, you'd need to make sure that you end
	// up with only one jump between the RNA and the protein

	// We want a fold tree that puts all the consecutive protein residues in one edge
	// and all consecutive RNA residues in another edge, with a jump between the
	// RNA and protein edge
	// Right now this will fail if there is e.g. RNA -> protein -> RNA or
	// protein -> RNA -> protein

	// We also need to make sure that linear_chainbreak is getting computed....

	kinematics::FoldTree ft;
	bool prev_RNA = false;
	bool prev_protein = false;
	vector1< core::Size > rna_protein_jumps;

	for ( core::Size i=1; i <= pose.size(); ++i ) {
		if ( pose.residue( i ).is_RNA() && prev_protein ) {
			rna_protein_jumps.push_back( i );
		} else if ( pose.residue( i ).is_protein() && prev_RNA ) {
			rna_protein_jumps.push_back( i );
		}
		if ( pose.residue( i ).is_RNA() ) {
			prev_RNA = true;
			prev_protein = false;
		} else if ( pose.residue( i ).is_protein() ) {
			prev_protein = true;
			prev_RNA = false;
		}
	}

	// Make the fold tree
	for ( core::Size i=1; i <= rna_protein_jumps.size(); ++i ) {
		if ( i == 1 ) { // add the first edge
			ft.add_edge(1, rna_protein_jumps[i] - 1, kinematics::Edge::PEPTIDE );
		}
		// add the jump
		ft.add_edge( rna_protein_jumps[i] - 1, rna_protein_jumps[i], i );
		if ( i == rna_protein_jumps.size() ) { // add the last edge
			ft.add_edge( rna_protein_jumps[i], pose.total_residue(), kinematics::Edge::PEPTIDE );
		}
	}

	return ft;

}



} //movers
} //denovo
} //rna
} //protocols
