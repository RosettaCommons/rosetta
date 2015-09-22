// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/motifs/IRCollection.cc
/// @brief Implmentation of interaction motifs

#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_Backrub.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/motifs/IRCollection.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/motif_utils.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <utility/string_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <iostream>
#include <map>
#include <algorithm>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

static THREAD_LOCAL basic::Tracer irt( "protocols.motifs.IRCollection", basic::t_info );

IRCollection::IRCollection() {}

IRCollection::IRCollection( core::pose::Pose & pose, MotifLibrary & motifs, utility::vector1< core::Size > const & build_sites )
{
	core::Size build_pos( 2 );

	// Loop over all build sites
	for ( core::Size target_index = 1, end_index = build_sites.size() ; target_index <= end_index ; ++target_index ) {
		core::Size target_pos( build_sites[ target_index ] );

		// Loop over motifs in libray, look for those that apply
		for ( MotifCOPs::const_iterator motif_itr_itr = motifs.begin(), end_itr = motifs.end() ; motif_itr_itr != end_itr ; ++motif_itr_itr ) {
			MotifCOP motif_itr( *motif_itr_itr );

			// Check for applicability
			if ( !motif_itr->apply_check( pose, target_pos ) ) continue;

			// irt << "Found a motif that applies!" << std::endl;
			// Build the set of inverse rotamers - note that I don't have a great idea for which phi/psi angles to
			// use for building.  The rotamer libraries are bb-dependent, so...
			bool is_it_forward( false );
			core::pack::rotamer_set::RotamerSetOP rot_set = motif_itr->build_inverted_rotamers( pose, target_pos, is_it_forward, build_pos );

			// Store the information
			motif_source_.push_back( motif_itr );
			motif_forward_.push_back( is_it_forward );
			target_positions_.push_back( target_pos );
			rotamer_sets_.push_back( rot_set );

		}
	}
}

core::Size
IRCollection::nirotamers() const
{
	core::Size accumulated_nirs( 0 );
	for ( core::Size index = 1, end_index = rotamer_sets_.size() ; index <= end_index ; ++index ) {
		accumulated_nirs += rotamer_sets_[ index ]->num_rotamers();
	}
	return accumulated_nirs;
}

void IRCollection::find_closest_backbone(
	core::pose::Pose & pose,
	protocols::loops::LoopsOP const flexible_regions,
	utility::vector1< core::Size > & closest_pos,
	utility::vector1< core::Real > & closest_rmsd
)
{
	core::Size this_nir( 0 );

	// Loop over rotamer sets in the collection
	for ( core::Size index = 1, end_index = rotamer_sets_.size() ; index <= end_index ; ++index ) {
		core::pack::rotamer_set::RotamerSetOP rotset( rotamer_sets_[ index ] );
		// Loop over rotamers in set
		for ( core::Size rindex = 1, end_rindex = rotset->num_rotamers() ; rindex <= end_rindex ; ++rindex ) {
			core::conformation::Residue const & this_rotamer( *rotset->rotamer( rindex ) );
			// Loop over backbone positions to consider

			this_nir++;

			for ( core::Size test_pos = 1, end_pos = pose.total_residue() ; test_pos <= end_pos ; ++test_pos ) {

				if ( !flexible_regions->is_loop_residue( test_pos ) ) {
					continue;
				}

				core::Real test_rmsd = backbone_stub_match( pose.residue( test_pos ), this_rotamer );
				if ( test_rmsd < closest_rmsd[ this_nir ] ) {
					closest_pos[ this_nir ] = test_pos;
					closest_rmsd[ this_nir ] = test_rmsd;
				}

			}
		}
	}
}

void IRCollection::incorporate_motifs( core::pose::Pose & pose, protocols::loops::LoopsOP const flexible_regions, utility::vector1< core::Size > & trim_positions )
{
	// Set up the initial state stuff and fire off the actual
	// recursive routine

	reset_unique_id();

	core::Size start_depth( 1 ); // Current recursion depth
	std::map< core::Size, MotifCOP > setpos; // previously incorporated positions
	std::map< core::Size, core::conformation::ResidueCOP > setpos_ir; // previously incorporated positions
	std::map< core::Size, bool > setpos_forward_info; // whether previous motifs are forward or reverse

	// makes loop residues alanine, except pro and gly, which stay the same
	mutate_loops_for_search( pose, *flexible_regions );
	mutate_position_vector_for_search( pose, trim_positions );

	// Relax away any clashes in the initial loop
	//core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep_design" ) );;
	core::scoring::ScoreFunctionOP score_fxn( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );
	core::scoring::methods::EnergyMethodOptions options( score_fxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	core::scoring::hbonds::HBondOptionsOP hb_options( new core::scoring::hbonds::HBondOptions );
	hb_options->exclude_DNA_DNA( false );
	hb_options->use_hb_env_dep( false );
	options.hbond_options( *hb_options );
	score_fxn->set_energy_method_options( options );
	protocols::loops::loop_mover::refine::LoopMover_Refine_Backrub pre_loop_refine( flexible_regions, score_fxn );
	//protocols::loops::loop_mover::refine::LoopMover_Refine_KIC pre_loop_refine( flexible_regions, score_fxn );
	pre_loop_refine.apply( pose );

	pose.dump_pdb( "alanine.pdb" );

	try_for_more( pose, flexible_regions, setpos, setpos_ir, setpos_forward_info, start_depth );

	return;
}

// This is the handler that tries to incorporate motifs into a flexible backbone region ad nauseum.
// It does the logic for which motifs to incorporate and spawns further branchings of the tree search
// recursively.  The one notable thing it does not do is the attempt at loop closure/relaxtion onto
// the motif itself.
void IRCollection::try_for_more(
	core::pose::Pose & start_pose,
	protocols::loops::LoopsOP const flexible_regions,
	std::map< core::Size, MotifCOP > setpos, // Note:  Passed by value (not by reference) intentionally
	std::map< core::Size, core::conformation::ResidueCOP > setpos_ir, // Note:  Passed by value (not by reference) intentionally
	std::map< core::Size, bool > setpos_forward_info, // Note:  Passed by value (not by reference) intentionally
	core::Size depth
)
{
	// Here we set up a number of parameters, cutoffs, etc. that should be data members of this
	// class.  The class itself should be renamed from IRCollection to MotifSearch or something
	// to denote that it does the motif incorporation and is not just a bag of inverse rotamers.

	core::Size const nirs( nirotamers() );

	utility::vector1< core::Size > closest_pos( nirs, 0 );
	utility::vector1< core::Real > closest_rmsd( nirs, 9999.0 );

	core::Real motif_rmsd_cutoff( basic::options::option[ basic::options::OptionKeys::motifs::close_enough ]() );
	core::Size max_depth( basic::options::option[ basic::options::OptionKeys::motifs::max_depth ]() );

	core::Size this_nir( 0 ); // IRs are stored in rotamer sets - this is maintained to count up in a flat array

	irt << "Entering search at recursion depth " << depth << std::endl;

	if ( depth <= max_depth ) {

		// Print out the positions with motifs already incorporated
		for ( std::map< core::Size, MotifCOP >::const_iterator itr = setpos.begin(), e_itr = setpos.end() ; itr != e_itr ; ++itr ) {
			irt << "Position " << itr->first << " is set." << std::endl;
		}

		// Find the closest backbone stubs for each inverse rotamer
		find_closest_backbone( start_pose, flexible_regions, closest_pos, closest_rmsd );

		// Loop through the full set of inverse rotamers
		for ( core::Size index = 1, end_index = rotamer_sets_.size() ; index <= end_index ; ++index ) {
			core::pack::rotamer_set::RotamerSetOP rotset( rotamer_sets_[ index ] );


			// Loop over rotamers in set
			for ( core::Size rindex = 1, end_rindex = rotset->num_rotamers() ; rindex <= end_rindex ; ++rindex ) {
				core::conformation::Residue const & this_rotamer( *rotset->rotamer( rindex ) );
				// Loop over backbone positions to consider

				this_nir++;
				core::Size const this_pos( closest_pos[ this_nir ] );

				// Here we disqualify positions that are already set, and those not in the defined flexible region
				if ( setpos.find( this_pos ) != setpos.end() ||
						!flexible_regions->is_loop_residue( this_pos ) ) {
					continue;
				}

				// Here we disqualify inverse rotamers that are too far away
				if ( closest_rmsd[ this_nir ] > motif_rmsd_cutoff ) {
					continue;
				}

				// Add a test to make sure this inverse rotamer doesn't clash with any of those previouisly incorporated


				// Attempt the incorporation for this rotamer
				irt << "Attempting to incorporate a motif at position " << this_pos << " with initial rmsd of " <<
					closest_rmsd[ this_nir ] << std::endl;

				core::pose::Pose pose;
				pose = start_pose;

				// Loop closure attempt
				if ( successful_loop_closure( pose, flexible_regions, setpos, setpos_ir, setpos_forward_info, this_rotamer, this_pos, motif_source_[ index ], motif_forward_[ index ]  ) ) {

					core::pose::Pose new_start_pose;
					new_start_pose = pose;

					std::map< core::Size, MotifCOP > new_setpos( setpos );
					new_setpos[ this_pos ] = motif_source_[ index ];
					std::map< core::Size, core::conformation::ResidueCOP > new_setpos_ir( setpos_ir );
					new_setpos_ir[ this_pos ] = this_rotamer.get_self_ptr();
					std::map< core::Size, bool > new_setpos_forward_info( setpos_forward_info );
					new_setpos_forward_info[ this_pos ] = motif_forward_[ index ];

					try_for_more( new_start_pose, flexible_regions, new_setpos, new_setpos_ir, new_setpos_forward_info, depth+1 );

				}
			}
		}
	}

	// Manipulate the pdb info to hide info about the motif in the pdb file

	// Score and output the file
	std::string unique_name( basic::options::option[ basic::options::OptionKeys::motifs::file_prefix ] );

	//  std::string unique_name( make_loop_pdb_name( setpos, unique_id() ) );
	increment_unique_id();
	std::string motif_filename( make_motif_filename( setpos, setpos_forward_info, start_pose ) );
	//std::string my_output_path( basic::options::option[ basic::options::OptionKeys::out::path::pdb ]() );
	std::string filename( unique_name + "_" + motif_filename + "_" + utility::to_string( unique_id() ) + ".pdb" );
	//std::string try_filename( my_output_path + "/" + unique_name + "_" + motif_filename + "_" + utility::to_string( unique_id() ) + ".pdb" );

	//core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep_design" ) );;
	core::scoring::ScoreFunctionOP score_fxn( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );
	core::scoring::methods::EnergyMethodOptions options( score_fxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	core::scoring::hbonds::HBondOptionsOP hb_options( new core::scoring::hbonds::HBondOptions );
	hb_options->exclude_DNA_DNA( false );
	hb_options->use_hb_env_dep( false );
	options.hbond_options( *hb_options );
	score_fxn->set_energy_method_options( options );

	start_pose.dump_scored_pdb( filename, *score_fxn );

	irt << "Finished search at recursion depth " << depth << std::endl;

	return;
}


bool
IRCollection::successful_loop_closure(
	core::pose::Pose & pose,
	protocols::loops::LoopsOP flexible_regions,
	std::map< core::Size, MotifCOP > & setpos,
	std::map< core::Size, core::conformation::ResidueCOP > & setpos_ir,
	std::map< core::Size, bool > & setpos_forward_info,
	core::conformation::Residue const & this_rotamer,
	core::Size const this_pos,
	MotifCOP this_motif,
	bool const this_forward_info
)
{
	using namespace core::scoring;

	core::Real coordinate_cst_weight( 10.0 );
	core::Real num_inverse_rotamers( 1.0 + setpos.size() ); // not an int because I need it for math

	//ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( "soft_rep_design" ) );;
	ScoreFunctionOP score_fxn( get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );
	methods::EnergyMethodOptions options( score_fxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	hbonds::HBondOptionsOP hb_options( new core::scoring::hbonds::HBondOptions );
	hb_options->exclude_DNA_DNA( false );
	hb_options->use_hb_env_dep( false );
	options.hbond_options( *hb_options );
	score_fxn->set_energy_method_options( options );
	score_fxn->set_weight( coordinate_constraint, coordinate_cst_weight );

	// Get bb constraints to this inverse rotamer

	constraints::ConstraintSetOP bb_cst_set( new constraints::ConstraintSet() );
	add_motif_bb_constraints( bb_cst_set, pose, this_pos, this_rotamer );

	// For previously set rotamers, set the sidechain constraints only
	for ( std::map< core::Size, MotifCOP >::const_iterator itr = setpos.begin(), e_itr = setpos.end() ; itr != e_itr ; ++itr ) {
		core::Size const prev_pos( itr->first );
		add_motif_sc_constraints( bb_cst_set, pose, prev_pos, *(setpos_ir[ prev_pos ]), itr->second, setpos_forward_info[ prev_pos ] );
	}

	pose.constraint_set( bb_cst_set );

	(*score_fxn)(pose);
	irt << "Before backbone refinement constraints score is " << pose.energies().total_energies()[ coordinate_constraint ] / num_inverse_rotamers << std::endl;

	// Perform some kind of loop relaxation
	protocols::loops::loop_mover::refine::LoopMover_Refine_Backrub loop_refine( flexible_regions, score_fxn );
	//protocols::loops::loop_mover::refine::LoopMover_Refine_KIC loop_refine( flexible_regions, score_fxn );

	loop_refine.apply( pose );

	// Evaluate the constraints energy to see how we did
	(*score_fxn)(pose);
	core::Real bb_constraint_check( pose.energies().total_energies()[ coordinate_constraint ] / num_inverse_rotamers );
	irt << "After backbone refinement constraints score is " << bb_constraint_check << std::endl;

	pose.remove_constraints();

	if ( bb_constraint_check > 1.0 ) {
		return false;
	}

	// Move the inverse rotamer onto the backbone
	pose.replace_residue( this_pos, this_rotamer, true );

	// Now put constraints on the actual atoms in the motif

	constraints::ConstraintSetOP sc_cst_set( new constraints::ConstraintSet() );
	add_motif_sc_constraints( sc_cst_set, pose, this_pos, this_rotamer, this_motif, this_forward_info );

	// For previously set rotamers, set their sidechain constraints as well
	for ( std::map< core::Size, MotifCOP >::const_iterator itr = setpos.begin(), e_itr = setpos.end() ; itr != e_itr ; ++itr ) {
		core::Size const prev_pos( itr->first );
		add_motif_sc_constraints( sc_cst_set, pose, prev_pos, *(setpos_ir[ prev_pos ]), itr->second, setpos_forward_info[ prev_pos ] );
	}

	pose.constraint_set( sc_cst_set );

	(*score_fxn)(pose);
	irt << "Before sidechain refinement constraints score is " << pose.energies().total_energies()[ coordinate_constraint ] / num_inverse_rotamers << std::endl;

	protocols::loops::loop_mover::refine::LoopMover_Refine_Backrub sc_loop_refine( flexible_regions, score_fxn );
	//protocols::loops::loop_mover::refine::LoopMover_Refine_KIC sc_loop_refine( flexible_regions, score_fxn );

	sc_loop_refine.apply( pose );

	// Evaluate the constraints energy to see how we did
	(*score_fxn)(pose);
	core::Real sc_constraint_check( pose.energies().total_energies()[ coordinate_constraint ] / num_inverse_rotamers );
	irt << "After sidechain refinement constraints score is " << sc_constraint_check << std::endl;

	pose.remove_constraints();

	if ( sc_constraint_check > 1.0 ) {
		return false;
	}

	return true;
}

std::string
IRCollection::make_motif_filename(
	std::map< core::Size, MotifCOP > & setpos,
	std::map< core::Size, bool > & setpos_forward_info,
	core::pose::Pose & pose
)
{
	std::string accum_string;
	for ( std::map< core::Size, MotifCOP >::const_iterator itr = setpos.begin(), e_itr = setpos.end()  ; itr != e_itr ; ++itr ) {
		core::Size const this_pos( itr->first );
		if ( setpos_forward_info[ this_pos ] ) {
			accum_string += itr->second->restype_name2() + utility::to_string( pose.pdb_info()->number( this_pos ) );
		} else {
			accum_string += itr->second->restype_name1() + utility::to_string( pose.pdb_info()->number( this_pos ) );
		}
	}
	return accum_string;
}

}
}
