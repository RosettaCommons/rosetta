// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kkappel/mutate_and_score_RNP.cc
/// @brief Make mutations to protein, RNA, or RNA/protein complex, update structure, score, and compare to wt
/// @author Kalli Kappel kappel@stanford.edu


//Unit Headers
//Package Headers
//Project Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/chemical/ChemicalManager.hh>
//Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
//Numeric Headers
#include <numeric/random/random.hh>
//C++ Headers
#include <iostream>
#include <fstream>

#include <devel/init.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/rigid/RigidBodyMover.hh>
//#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// New options for this application
using namespace basic::options::OptionKeys;

OPT_KEY( String, mutfile )
OPT_KEY( Boolean, move_backbone )
OPT_KEY( Boolean, move_protein_backbone )
OPT_KEY( Boolean, min_jumps )
OPT_KEY( String, out_prefix )
OPT_KEY( String, unbound_protein )
OPT_KEY( String, unbound_RNA )
OPT_KEY( FileVector, unbound_wtRNA_ens )
OPT_KEY( Boolean, mutate_unbound_RNA )
OPT_KEY( Integer, prot_offset_num )
OPT_KEY( Integer, RNA_offset_num )
OPT_KEY( Integer, protein_start_res )
OPT_KEY( Integer, protein_end_res )
OPT_KEY( Real, relax_cutoff_dist )
OPT_KEY( Boolean, min_only )
OPT_KEY( Boolean, dump_bound_rna )
OPT_KEY( String, bound_rna_dump_tag )
OPT_KEY( Integer, protein_pack_reps )
OPT_KEY( Integer, Nreps )

static basic::Tracer TR( "apps.public.rnp_ddg.rnp_ddg" );

/////////////////////////////////////////////////////////////////////////////////
// Get mutations from input file and store in private variable
std::pair< utility::vector1< utility::vector1< std::pair< core::Size, char > > >, utility::vector1< std::string > >
get_mutations_from_file(
	std::string mut_file,
	core::pose::Pose const & pose
) {
	utility::vector1< utility::vector1< std::pair< core::Size, char > > > list_of_mutations;
	utility::vector1< std::string > list_of_input_files;

	// Read from file
	std::ifstream mut_input;
	mut_input.open( mut_file.c_str() );

	if ( !mut_input.is_open() ) {
		TR << "ERROR: Can't open the mutfile!" << std::endl;
		exit (1);
	}

	while ( !mut_input.eof() ) {
		std::string elem;
		mut_input >> elem;
		if ( elem == "total" ) {
			core::Size num_mutations, check_num_mutations, resnum;
			mut_input >> num_mutations;
			mut_input >> check_num_mutations;
			if ( num_mutations != check_num_mutations ) {
				TR << "ERROR: Please check your mutations file!! Are you sure about the total number of mutations?" << std::endl;
				exit (1);
			}
			utility::vector1< std::pair< core::Size, char > > mutations;
			for ( core::Size i = 1; i <= num_mutations; ++i ) {
				char wt_res, mut_res;
				mut_input >> wt_res;
				mut_input >> resnum;
				mut_input >> mut_res;
				resnum += basic::options::option[ prot_offset_num ](); //default = 0 //TODO: this assumes the RNA mutation!
				resnum += basic::options::option[ RNA_offset_num ](); //default = 0
				// Check that the wildtype residue is correct
				TR << "Trying to mutate " << wt_res << resnum << " to " << mut_res << std::endl;
				TR << "You say wt is " << wt_res << " and it is really " << pose.residue( resnum ).name1() << std::endl;
				runtime_assert( pose.residue( resnum ).name1() == wt_res );
				std::pair< core::Size, char > mutation;
				mutation = std::make_pair( resnum, mut_res );
				mutations.push_back( mutation );
			}
			list_of_mutations.push_back( mutations );
		} else if ( elem == "input_file" ) { // elem was not equal to "total"
			// if input_file was specified, don't make mutations, make note of the file name
			std::string file_name;
			mut_input >> file_name;
			// TODO: should check that this contains ".pdb"
			list_of_input_files.push_back( file_name );
		}
	}

	std::pair< utility::vector1< utility::vector1< std::pair< core::Size, char > > >, utility::vector1< std::string > > all_input;
	all_input = std::make_pair( list_of_mutations, list_of_input_files );

	// Check that if any input files were specified, they were specified for all mutants
	if ( list_of_input_files.size() > 0 and (list_of_mutations.size() != list_of_input_files.size()) ) {
		TR << "ERROR: you didn't specify enough input files in the mutfile, must specify all or none!" << std::endl;
		TR << "Use 'input_file none' for mutants without input structures" << std::endl;
	}

	return all_input;
}


/////////////////////////////////////////////////////////////////////////////////
//Find the residues around the mutation
//Return a vector of vectors of residues near the mutations?
utility::vector1< core::Size > find_residues_around_mutation(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > mutant_residues )
{
	using namespace basic::options;
	utility::vector1< core::Size > neighbors;
	// Totally arbitary, set the cutoff distance for neighbors to be 10A
	core::Real const cutoff = option[ relax_cutoff_dist ]();

	// Loop through all the residues in the pose to find those that are near
	// the mutations
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		// If this residue is one of the mutation positions just add it
		// to the list and continue
		for ( core::Size j = 1; j <= mutant_residues.size(); ++j ) {
			if ( i == mutant_residues[ j ] ) {  // The residue is in the list of mutations
				// Check that this residues isn't already in the list of neighbors
				if ( !neighbors.has_value( i ) ) {
					neighbors.push_back( i );
				}
				continue;
			}
		}
		core::Vector pos1, pos2;
		// Now if the residue is not one of the mutation positions
		if ( pose.residue( i ).is_protein() ) {
			pos1 = pose.residue( i ).xyz( " CA " );
			//Vector pos1( pose.residue( i ).xyz( " CA " ) );
		} else if ( pose.residue( i ).is_RNA() ) {
			pos1 = pose.residue( i ).xyz( " P  " );
		} else { // if it's not RNA or protein
			TR << "WARNING: Residues other than RNA or protein are not going to be minimized!!" << std::endl;
			continue;
		}
		// Loop through the mutations and calculate distances
		for ( core::Size j = 1; j <= mutant_residues.size(); ++j ) {
			if ( pose.residue( mutant_residues[ j ] ).is_protein() ) {
				pos2 = pose.residue( mutant_residues[ j ] ).xyz( " CA " );
			} else if ( pose.residue( mutant_residues[ j ] ).is_RNA() ) {
				pos2 = pose.residue( mutant_residues[ j ] ).xyz( " P  " );
			} else { //the mutation position is not RNA or protein!
				TR << "ERROR: Cannot handle mutating something other than RNA or protein!!" << std::endl;
				exit( 1 );
			}
			core::Real distance = (pos1 - pos2).length();
			if ( distance < cutoff ) {
				// Check that this residues isn't already in the list of neighbors
				if ( !neighbors.has_value( i ) ) {
					neighbors.push_back( i );
				}
			}
		}
	}

	return neighbors;
}

/////////////////////////////////////////////////////////////////////////////////
// This version just takes a list of mutated positions, rather than a list of mutations
utility::vector1< core::Size > find_residues_around_mutation(
	core::pose::Pose const & pose,
	utility::vector1< std::pair< core::Size, char > > const mutations )
{
	utility::vector1< core::Size > mutant_residues;
	for ( core::Size i = 1; i<=mutations.size(); ++i ) {
		mutant_residues.push_back(mutations[i].first);
	}
	return find_residues_around_mutation( pose, mutant_residues );
}

/////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP setup_score_function()
{
	core::scoring::ScoreFunctionOP sfxn;
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		// Get the score function specified by the user on the command line
		sfxn = core::scoring::get_score_function();
	} else { // Set a default score function
		sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "flicker2015_sol2_fa_elec3.wts" );
	}

	return sfxn;
}

/////////////////////////////////////////////////////////////////////////////////
void simple_minimization(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP sfxn,
	utility::vector1< core::Size > const neighbor_residues,
	bool exclude_separation_jump = false ) // could instead provide a list of jumps to exclude
{
	using namespace core::scoring;
	// Save initial atom_pair_constraint weight
	core::Real init_cst_wt = sfxn->get_weight( atom_pair_constraint );

	// Set atom_pair constraint weight to 1
	sfxn->set_weight( atom_pair_constraint, 1.0 );

	core::optimization::MinimizerOptions min_opts( "dfpmin_armijo_nonmonotone", 0.000025, true );
	// Used in stepwise, so I'll use it here as well
	// This slows things down considerably
	min_opts.nblist_auto_update( true );
	core::optimization::AtomTreeMinimizer atm;
	core::kinematics::MoveMap mm;
	// Loop through the list of neighbor residues and set the movemap to be true for them
	for ( core::Size i = 1; i <= neighbor_residues.size(); ++i ) {
		if ( basic::options::option[ move_backbone ]() ) {
			if ( pose.residue( neighbor_residues[i] ).is_protein() && basic::options::option[ move_protein_backbone ]() ) {
				mm.set_bb( neighbor_residues[i], true );
			} else if ( pose.residue( i ).is_RNA() ) {
				mm.set_bb( neighbor_residues[i], true );
			}
		}
		mm.set_chi( neighbor_residues[i], true );
		mm.set_jump( true );
		if ( exclude_separation_jump ) {
			//TODO: Assuming that jump 1 is the separation jump -- not great
			mm.set_jump( 1 /* separation jump */, false );
		}
	}
	// Option to exclude minimization of jumps
	if ( !basic::options::option[ min_jumps ]() ) {
		mm.set_jump( false );
	}
	TR << "Doing minimization..." << std::endl;
	atm.run( pose, mm, *sfxn, min_opts );
	TR << "Finished minimization!" << std::endl;

	// Reset the atom_pair_constraint weight
	sfxn->set_weight( atom_pair_constraint, init_cst_wt );
}

/////////////////////////////////////////////////////////////////////////////////
// Overload this function so that it can take a list of neighbor residues or not
void simple_minimization(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP sfxn,
	bool exclude_separation_jump = false ) // could instead provide a list of jumps to exclude
{
	utility::vector1< core::Size > neighbor_residues;
	// All the residues in the pose will be added to the list of neighbor residues
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		neighbor_residues.push_back( i );
	}
	simple_minimization( pose, sfxn, neighbor_residues, exclude_separation_jump );
}

/////////////////////////////////////////////////////////////////////////////////
void unbound_protein_minimization(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP sfxn )
{
	using namespace core::scoring;
	// Save initial atom_pair_constraint weight
	core::Real init_cst_wt = sfxn->get_weight( atom_pair_constraint );

	// Set atom_pair constraint weight to 1
	sfxn->set_weight( atom_pair_constraint, 1.0 );

	core::optimization::MinimizerOptions min_opts( "dfpmin_armijo_nonmonotone", 0.000025, true );
	// Used in stepwise, so I'll use it here as well
	// This slows things down considerably
	min_opts.nblist_auto_update( true );
	core::optimization::AtomTreeMinimizer atm;
	core::kinematics::MoveMap mm;
	// Loop through the list of neighbor residues and set the movemap to be true for them
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		// Do not move backbone at all for now
		mm.set_chi( i, true );
		mm.set_jump( true ); // should only be jumps between protein residues
	}
	// Option to exclude minimization of jumps
	if ( !basic::options::option[ min_jumps ]() ) {
		mm.set_jump( false );
	}
	TR << "Doing minimization..." << std::endl;
	atm.run( pose, mm, *sfxn, min_opts );
	TR << "Finished minimization!" << std::endl;

	// Reset the atom_pair_constraint weight
	sfxn->set_weight( atom_pair_constraint, init_cst_wt );
}

/////////////////////////////////////////////////////////////////////////////////

//Check out setup_pack_task in StepWisePacker for example packing RNA
void simple_packing(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP sfxn,
	utility::vector1< core::Size > const neighbor_residues )
{
	using namespace core::scoring;
	core::pack::task::PackerTaskOP pack_task( core::pack::task::TaskFactory::create_packer_task( pose ));
	pack_task->restrict_to_repacking();
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		// Only pack the neighbor residues
		if ( neighbor_residues.has_value( i ) ) {
			pack_task->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
			pack_task->nonconst_residue_task(i).or_include_current( true );
			if ( pose.residue( i ).is_protein() ) {
				pack_task->nonconst_residue_task( i ).or_ex1( true );
				pack_task->nonconst_residue_task( i ).or_ex2( false );
			} else if ( pose.residue( i ).is_RNA() ) {
				// Need or_ex1 true to get more rotamers than just current
				pack_task->nonconst_residue_task( i ).or_ex1( true );
				pack_task->nonconst_residue_task( i ).or_ex4( true );
			}
		} else {
			// If not in the list of neighbors, don't pack
			pack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}
	TR << "Doing packing..." << std::endl;
	core::pack::pack_rotamers( pose, *sfxn, pack_task );
	TR << "Finished packing!" << std::endl;
}
/////////////////////////////////////////////////////////////////////////////////
void simple_packing(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP sfxn )
{
	utility::vector1< core::Size > neighbor_residues;
	// All the residues in the pose will be added to the list of neighbor residues
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		neighbor_residues.push_back( i );
	}
	simple_packing( pose, sfxn, neighbor_residues );
}

/////////////////////////////////////////////////////////////////////////////////
// Want to do the mutations differently for RNA vs protein
// Protein use the packer, RNA use mutate position
void make_mutations(
	core::pose::Pose & pose_to_mutate,
	core::scoring::ScoreFunctionOP sfxn,
	utility::vector1< std::pair< core::Size, char > > const mutations,
	utility::vector1< core::Size > const neighbor_residues,
	bool override_did_mutate_and_pack )
{
	bool did_mutate = false;
	utility::vector1< core::Size > protein_mut_positions;
	utility::vector1< char > protein_mut_types;
	for ( core::Size i = 1; i <= mutations.size(); ++i ) {
		if ( pose_to_mutate.residue( mutations[i].first ).is_RNA() ) {
			core::pose::rna::mutate_position( pose_to_mutate, mutations[i].first, mutations[i].second );
			did_mutate = true;
		} else if ( pose_to_mutate.residue( mutations[i].first ).is_protein() ) {
			protein_mut_positions.push_back( mutations[i].first );
			protein_mut_types.push_back( mutations[i].second );
		} else {
			// This mutation is not RNA or protein!
			TR << "WARNING!!! Cannot handle a mutation of a residue type other than RNA or protein!" << std::endl;
			exit (1);
		}
	}
	// Do the protein mutations, if any
	if ( protein_mut_positions.size() > 0 ) {
		//Create a pack task for the protein mutations
		core::pack::task::PackerTaskOP pack_task( core::pack::task::TaskFactory::create_packer_task( pose_to_mutate ));
		// Loop through all the residues in the pose
		for ( core::Size r = 1; r <= pose_to_mutate.total_residue(); ++r ) {
			core::Size index = std::find( protein_mut_positions.begin(), protein_mut_positions.end(), r ) - protein_mut_positions.begin();
			if ( index < protein_mut_positions.size() ) { // this is a protein mutation position
				core::chemical::AA mut_type( core::chemical::aa_from_oneletter_code( protein_mut_types[ index+1 ] ));
				utility::vector1<bool> restrict_to_aa( 20, false );
				restrict_to_aa[ (core::Size)mut_type ] = true;
				TR.Debug << restrict_to_aa << std::endl;
				pack_task->nonconst_residue_task( r ).restrict_absent_canonical_aas( restrict_to_aa );
				pack_task->nonconst_residue_task( r ).or_include_current( true );
				pack_task->nonconst_residue_task( r ).or_ex1( true );
				pack_task->nonconst_residue_task( r ).or_ex2( true );
				pack_task->nonconst_residue_task(r).print_allowed_types( TR.Debug );
			} else {
				pack_task->nonconst_residue_task( r ).restrict_to_repacking();
				pack_task->nonconst_residue_task( r ).prevent_repacking();
			}
		}
		// Do the packing (i.e. make the mutations)
		core::pack::pack_rotamers( pose_to_mutate, *sfxn, pack_task );
		did_mutate = true;
	}

	if ( did_mutate || override_did_mutate_and_pack ) {
		// Minimizing immediately sometimes creates very bad structures! 04-22-16
		// But sometimes packing immediately after mutating creates weird structures
		// So take the minimum scoring structure of (min, pack, min) and (pack, min)

		// Make a hard copy of the pose
		core::pose::Pose mutant_pose_1 = pose_to_mutate;
		// Just pack and then min the copy pose:
		if ( !basic::options::option[ min_only ]() ) {
			simple_packing( mutant_pose_1, sfxn, neighbor_residues );
		}
		simple_minimization( mutant_pose_1, sfxn, neighbor_residues );

		// Min, pack, then min the original mutant pose
		simple_minimization( pose_to_mutate, sfxn, neighbor_residues );
		if ( !basic::options::option[ min_only ]() ) {
			simple_packing( pose_to_mutate, sfxn, neighbor_residues );
		}
		simple_minimization( pose_to_mutate, sfxn, neighbor_residues );

		// Figure out which one scores better
		core::Real score_min_pack_min = (*sfxn)( pose_to_mutate );
		core::Real score_pack_min = (*sfxn)( mutant_pose_1 );
		// If pack_min scores better, copy it to pose_to_mutate (which will exist outside this function)
		if ( score_pack_min < score_min_pack_min ) {
			TR.Debug << "score_pack_min is lower" << std::endl;
			pose_to_mutate = mutant_pose_1;
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////
core::kinematics::FoldTree get_simple_RNA_protein_fold_tree(
	core::pose::Pose const & pose )
{

	utility::vector1< core::Size > protein_start_residues, RNA_start_residues, all_start_residues;
	core::kinematics::FoldTree separation_foldtree;

	// Loop through pose and find residues that start each block of protein or RNA
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( i == 1 && pose.residue( i ).is_protein() ) {
			protein_start_residues.push_back( i );
			all_start_residues.push_back( i );
		} else if ( i ==1 && pose.residue( i ).is_RNA() ) {
			RNA_start_residues.push_back( i );
			all_start_residues.push_back( i );
		} else if ( pose.residue( i ).is_RNA() && pose.residue( i-1 ).is_protein() ) {
			RNA_start_residues.push_back( i );
			all_start_residues.push_back( i );
		} else if ( pose.residue( i ).is_protein() && pose.residue( i-1).is_RNA() ) {
			protein_start_residues.push_back( i );
			all_start_residues.push_back( i );
		}
	}

	// Now set up the new fold tree
	// Add all the edges using the list of all the start residues
	for ( core::Size i = 1; i < all_start_residues.size(); ++i ) {
		separation_foldtree.add_edge( all_start_residues[i], all_start_residues[i+1]-1,
			core::kinematics::Edge::PEPTIDE ); // I think we just use peptide even though it's sometime RNA?
	}
	// Then add the last edge
	separation_foldtree.add_edge( all_start_residues.back(), pose.total_residue(), core::kinematics::Edge::PEPTIDE );

	// Now add the jumps; these should go from the first residue to the second start residue
	// And from the first residue to all other start residues of the same type
	// And from the second start residue to all other start residues of the same type

	// Add the jump from the first residue to the second start residue
	if ( protein_start_residues[ 1 ] == 1 ) {
		separation_foldtree.add_edge( protein_start_residues[ 1 ], RNA_start_residues[ 1 ], 1 /* jump # */ );
	} else {
		separation_foldtree.add_edge( RNA_start_residues[ 1 ], protein_start_residues[ 1 ], 1 /* jump # */ );
	}

	// Add the jumps between all protein residues
	for ( core::Size i = 1; i < protein_start_residues.size(); ++i ) {
		separation_foldtree.add_edge( protein_start_residues[ i ], protein_start_residues[ i+1 ], i+1 /* jump # */ );
	}

	// Add the jumps between all RNA residues
	for ( core::Size i = 1; i < RNA_start_residues.size(); ++i ) {
		separation_foldtree.add_edge( RNA_start_residues[ i ], RNA_start_residues[ i+1 ], i+1 /* jump # */ );
	}

	TR.Debug << "The reordered separation fold tree" << std::endl;
	separation_foldtree.reorder( 1 );
	separation_foldtree.show( TR.Debug );

	return separation_foldtree;
}

/////////////////////////////////////////////////////////////////////////////////
void delete_protein_from_pose( core::pose::Pose & pose ) {

	// Delete protein from pose

	utility::vector1< std::pair< core::Size, core::Size > > protein_start_and_end_residues;
	core::Size start = 0;
	core::Size end = 0;
	// Figure out the protein start and end residues
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( pose.residue( i ).is_protein() && start == 0 ) {
			start = i;
		} else if ( !pose.residue( i ).is_protein() ) {
			if ( start != 0 ) {
				end = i-1;
				protein_start_and_end_residues.push_back( std::make_pair( start, end ) );
				start = 0;
			}
		}
		if ( (i == pose.total_residue()) && (start > end) ) {
			end = i;
			protein_start_and_end_residues.push_back( std::make_pair( start, end ) );
		}
	}

	// Now delete all the protein residues from the pose
	core::Size offset = 0;
	for ( core::Size i=1; i<=protein_start_and_end_residues.size(); ++i ) {
		pose.delete_residue_range_slow( protein_start_and_end_residues[i].first - offset, protein_start_and_end_residues[i].second - offset );
		offset += protein_start_and_end_residues[i].second - protein_start_and_end_residues[i].first + 1;
	}

}
/////////////////////////////////////////////////////////////////////////////////

void delete_RNA_from_pose( core::pose::Pose & pose ) {

	// Delete RNA from pose

	utility::vector1< std::pair< core::Size, core::Size > > rna_start_and_end_residues;
	core::Size start = 0;
	core::Size end = 0;
	// Figure out the protein start and end residues
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( pose.residue( i ).is_RNA() && start == 0 ) {
			start = i;
		} else if ( !pose.residue( i ).is_RNA() ) {
			if ( start != 0 ) {
				end = i-1;
				rna_start_and_end_residues.push_back( std::make_pair( start, end ) );
				start = 0;
			}
		}
		if ( (i == pose.total_residue()) && (start > end) ) {
			end = i;
			rna_start_and_end_residues.push_back( std::make_pair( start, end ) );
		}
	}

	// Now delete all the protein residues from the pose
	core::Size offset = 0;
	for ( core::Size i=1; i<=rna_start_and_end_residues.size(); ++i ) {
		pose.delete_residue_range_slow( rna_start_and_end_residues[i].first - offset, rna_start_and_end_residues[i].second - offset );
		offset += rna_start_and_end_residues[i].second - rna_start_and_end_residues[i].first + 1;
	}

}

/////////////////////////////////////////////////////////////////////////////////
std::string make_mutation_tag(
	utility::vector1< std::pair< core::Size, char > > mutations )
{
	// Make a tag for the silent struct and score for the specific mutant
	std::string mutation_tag = basic::options::option[ out_prefix ]();
	if ( mutations.size() == 0 ) {
		mutation_tag.append( "wildtype" );
	}
	//mutation_tag = "wildtype"; }
	for ( core::Size i = 1; i <= mutations.size(); ++i ) {
		std::stringstream ss;
		ss << mutations[i].first;
		mutation_tag.append( ss.str() );
		mutation_tag.push_back( mutations[i].second );
	}

	return mutation_tag;
}

/////////////////////////////////////////////////////////////////////////////////
// Obviously this isn't the most accurate computation of binding energy
// ddG_bind = dG_Complex - dG_Receptor - dG_Ligand
core::Real calculate_binding_energy(
	core::scoring::ScoreFunctionOP & sfxn,
	core::pose::Pose & pose, // TODO: actually kept const because copying to unbound_pose...
	utility::vector1< std::pair< core::Size, char > > mutations )
{

	using namespace basic::options;
	using namespace core::scoring;
	using namespace core::import_pose::pose_stream;


	core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// TODO: Be sure that there aren't any constraints on!
	// Initialize scores to 0
	core::Real complex_score( 0.0);
	core::Real unbound_score( 0.0); // This may never get reset if protein_pack_reps=0

	core::Real ddG_bind( 0.0);
	// Turn off csts for scoring, TODO: might need to do this in a smarter way if these csts
	// become part of the normal score, or if we use something other than atom_pair_constraints
	core::Real initial_cst_wt = sfxn->get_weight( atom_pair_constraint );
	sfxn->set_weight( atom_pair_constraint, 0.0 );
	complex_score = (*sfxn)(pose);
	TR << "Score of the complex: " << complex_score << std::endl;

	core::pose::Pose unbound_pose = pose;
	core::kinematics::FoldTree separation_foldtree, saved_foldtree;
	saved_foldtree = unbound_pose.fold_tree();

	// Calculate the "interaction energy" of the complex and print this out
	core::Real interaction_E, unbound_E_norelax;

	separation_foldtree = get_simple_RNA_protein_fold_tree( unbound_pose );

	// Set the fold tree of the unbound pose to the new separation fold tree
	unbound_pose.fold_tree( separation_foldtree );

	// Now the first jump is the correct one to use for the separation of the RNA and protein
	protocols::rigid::RigidBodyTransMoverOP separate_complex( new protocols::rigid::RigidBodyTransMover( unbound_pose, 1 /* jump # */) );
	separate_complex->step_size( 1000.0 );
	separate_complex->apply( unbound_pose );

	// Just score the separated complex without any relaxation of structures
	unbound_E_norelax = (*sfxn)(unbound_pose);
	interaction_E = complex_score - unbound_E_norelax;

	TR << "Interaction energy of complex: " << interaction_E << std::endl;


	// Impose a temporary fold tree for the separation of the complex
	// Neccessary for example if you input PDB with RNA chain 1 + RNA chain 2 + protein
	// (for this you would need to use jump number 2 to separate RNA and protein)
	// This will assume that the complex you're using to calculate ddG is RNA/protein
	// Could give a second chance by trying to figure it out like this, and then if there's only RNA
	// or only protein, then have to assume chains are representative of the two partners and use jump = 1
	// for separation

	if ( option[ unbound_protein ].user() && option[ unbound_wtRNA_ens ].user() ) {
		// unbound_wtRNA_ens can include just one or many structures
		// This method does not work well (tested on MS2 initial set of RNA mutations)

		std::string protein_input = option[ unbound_protein ]();
		utility::vector1< std::string > RNA_input = option[ unbound_wtRNA_ens ]();
		PoseInputStreamOP protein_input_stream, rna_input_stream;
		core::pose::Pose protein_pose, RNA_pose;
		protein_input_stream = PoseInputStreamOP( new PDBPoseInputStream( protein_input ) );
		rna_input_stream = PoseInputStreamOP( new PDBPoseInputStream( RNA_input ) );
		protein_input_stream->fill_pose( protein_pose, *rsd_set );
		core::Real unbound_RNA_score, unbound_protein_score;
		// the ~partition function for a given secondary structure, assuming we are loading in all microstates
		core::Real partition_function_SS = 0.0;
		core::Real kT = 0.59; //Boltzmann constant at room temperature
		core::Size i = 0;
		while ( rna_input_stream->has_another_pose() ) {
			i += 1;
			rna_input_stream->fill_pose( RNA_pose, *rsd_set );
			// Do mutations of the unbound RNA, if specified, otherwise, just keep original seq
			if ( option[ mutate_unbound_RNA ]() ) {
				// Fix numbering of mutations (TODO: definitely a better way to do this)
				utility::vector1< std::pair< core::Size, char > > mut_rna_num;
				for ( core::Size n = 1; n <= mutations.size(); ++n ) {
					mut_rna_num.push_back( mutations[n] );
					mut_rna_num[n].first -= option[ prot_offset_num ]();
				}
				utility::vector1< core::Size > RNA_neighbor_residues =
					find_residues_around_mutation( RNA_pose, mut_rna_num );
				make_mutations( RNA_pose, sfxn, mut_rna_num, RNA_neighbor_residues, false );
				// Do another minimization of the full pose
				simple_minimization( RNA_pose, sfxn );
			}
			core::Real unbound_RNA_E = (*sfxn)( RNA_pose );
			TR << "Unbound RNA for ensemble member " << i << ": " << unbound_RNA_E << std::endl;
			partition_function_SS += std::exp( -1.0 * unbound_RNA_E / kT );
		}
		unbound_RNA_score = -1.0 * kT * std::log( partition_function_SS );
		unbound_protein_score = (*sfxn)( protein_pose );
		TR << "Unbound protein score: " << unbound_protein_score << std::endl;
		TR << "Total Unbound RNA score: " << unbound_RNA_score << std::endl;
		unbound_score = unbound_RNA_score + unbound_protein_score;
	} else if ( option[ unbound_protein ].user() && !option[ unbound_wtRNA_ens ].user() ) {
		// This method is only useful if you actually want to calculate ddG for only RNA mutations
		// by taking complex_score - protein_score(unbound) - RNA_score(unbound)
		// Then it will reduce some noise in the calculations by using just one unbound protein score/
		// structure for all mutants
		// Using -protein_pack_reps 0 is probably better though (it's not really necessary to calculate
		// the unbound protein score if it is assumed to be the same for all mutants) -- saves some computation

		core::Real unbound_RNA_score, unbound_protein_score;
		// Just read the protein pose and score for the unbound protein score
		core::pose::Pose protein_pose, RNA_pose;
		std::string protein_input = option[ unbound_protein ]();

		PoseInputStreamOP protein_input_stream;
		protein_input_stream = PoseInputStreamOP( new PDBPoseInputStream( protein_input ) );
		protein_input_stream->fill_pose( protein_pose, *rsd_set );

		unbound_protein_score = (*sfxn)( protein_pose );

		// Make a copy of the unbound pose so that we're not screwing anything up
		RNA_pose = unbound_pose;

		//  // Delete the protein residues from the pose and then minimize to get the unbound RNA score
		delete_protein_from_pose( RNA_pose );

		// If specified, dump out this RNA conf
		if ( option[ dump_bound_rna ]() ) {
			RNA_pose.dump_pdb( "RNA_bound_conf_" + option[ bound_rna_dump_tag ]() + ".pdb" );
		}
		simple_minimization( RNA_pose, sfxn );
		unbound_RNA_score = (*sfxn)( RNA_pose );

		unbound_score = unbound_RNA_score + unbound_protein_score;

		TR << "Unbound RNA score: " << unbound_RNA_score << std::endl;
		TR << "Unbound protein score: " << unbound_protein_score << std::endl;

	} else {
		// TODO: Want to do packing and minimization just of the interface residues?
		// For now just pack and minimize the whole structure
		// If want to calculate more accurate unbound protein score, then repeat this many times
		// and take the average score
		// If there aren't any protein mutations and don't actually care about the unbound protein
		// score then protein_pack_reps can be set to 0
		utility::vector1< core::Real > unbound_protein_scores;
		core::Real total_unbound_protein( 0.0 );
		core::Size pack_reps( option[ protein_pack_reps ]() );
		for ( core::Size i=1; i<= pack_reps; ++i ) { // Default 1 rep
			simple_minimization( unbound_pose, sfxn, true /*exclude separation jump*/ ); //TODO: Minimize the whole structure? I think we want to?
			if ( !option[min_only]() ) {
				simple_packing( unbound_pose, sfxn );
				simple_minimization( unbound_pose, sfxn, true /*exclude separation jump*/ ); //TODO: Minimize the whole structure? I think we want to?
			}

			// Score the unbound components
			unbound_score = (*sfxn)(unbound_pose);

			// Scores of individual components
			core::Real unbound_protein_score, unbound_RNA_score;
			core::pose::Pose RNA_pose = unbound_pose;

			// Delete the protein residues from the pose
			delete_protein_from_pose( RNA_pose );

			unbound_RNA_score = (*sfxn)( RNA_pose );
			unbound_protein_score = unbound_score - unbound_RNA_score;
			TR << "Unbound RNA score: " << unbound_RNA_score << std::endl;
			TR << "Unbound protein score: " << unbound_protein_score << std::endl;
			// append the protein score to the list
			unbound_protein_scores.push_back( unbound_protein_score );
			total_unbound_protein+=unbound_protein_score;
		}
		// Find the average unbound protein score
		if ( pack_reps > 1 ) {
			core::Real avg_unbound_protein_score( total_unbound_protein / pack_reps );
			TR << "Average unbound protein score: " << avg_unbound_protein_score << std::endl;
			// TODO: Find and print the standard deviation (?)
		}
	}
	TR << "Score of the unbound: " << unbound_score << std::endl;
	// Subtract unbound from complex to get ddG_bind
	ddG_bind = complex_score - unbound_score;
	TR << "ddG_bind is: " << ddG_bind << std::endl;

	// Output binding energy to some file

	// Create this ddG_silent_struct to be able to take advantage of the silent file
	// io code, this will neatly write out the differences in all the score terms

	// First update the full_model_info so the output sequence in the silent file is
	// actually correct
	core::pose::full_model_info::update_full_model_info_from_pose( pose );
	core::pose::full_model_info::update_full_model_info_from_pose( unbound_pose );

	// Make a tag for the silent struct and score for the specific mutant
	std::string mutation_tag = make_mutation_tag( mutations );

	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData ddG_silent_file_data( opts );
	// Use the bound structure to create the ddG_silent_struct, but this is arbitrary
	// because I will only use this for writing out scores, which are updated below
	// to be equal to the difference of the bound and unbound scores
	core::io::silent::BinarySilentStruct ddG_silent_struct( opts, pose, mutation_tag );
	std::string const ddG_silent_file = mutation_tag + "_ddG_bind.sc";

	//TODO: Need to fix this if input unbound rna and protein poses! (otherwise score file is 0's except total score)
	// Update the scores in the ddG_silent_struct to be equal to the difference of the
	// bound and unbound scores ( modification of code starting at line 651 in silent/SilentStruct.cc)
	core::scoring::EnergyMap const emap_bound = pose.energies().total_energies();
	core::scoring::EnergyMap const emap_unbound = unbound_pose.energies().total_energies();
	core::scoring::EnergyMap const wts  = pose.energies().weights();
	core::scoring::EnergyMap::const_iterator emap_unbound_iter, emap_bound_iter, wts_iter;
	// Iterate through the score terms for the bound and unbound structures
	for ( emap_bound_iter = emap_bound.begin(), wts_iter = wts.begin(), emap_unbound_iter = emap_unbound.begin();
			emap_bound_iter != emap_bound.end() && wts_iter != wts.end() && emap_unbound_iter != emap_unbound.end();
			++emap_bound_iter && ++wts_iter && ++emap_unbound_iter ) {
		if ( *wts_iter != 0.0 ) { // get terms w/ non-zero wts
			core::scoring::ScoreType sc_type
				= core::scoring::ScoreType( emap_bound_iter - emap_bound.begin() + 1 );
			std::string sc_name = core::scoring::name_from_score_type( sc_type );
			// Append the difference of the bound and the unbound scores
			ddG_silent_struct.add_energy( sc_name, *emap_bound_iter - *emap_unbound_iter, *wts_iter );
		}
	}
	// Set the total score term (= to binding energy)
	// Need to do this because we didn't iterate over this term above
	ddG_silent_struct.add_energy( "score", ddG_bind, 1.0 ); // this overwrites "score"
	// Write out the scores to a file, writes out ddG for each score term
	ddG_silent_file_data.write_silent_struct( ddG_silent_struct, ddG_silent_file, true /* score only */ );

	// Write a silent file with the bound and unbound structures
	std::string const mutant_struct_file = mutation_tag + "_structs.out"; //TODO: allow user to specify file name
	core::io::silent::SilentFileData mutant_silent_data( opts );
	core::io::silent::BinarySilentStruct mutant_silent_struct_bound( opts, pose, mutation_tag + "_bound" );
	core::io::silent::BinarySilentStruct mutant_silent_struct_unbound( opts, unbound_pose, mutation_tag + "_unbound" );
	mutant_silent_data.write_silent_struct( mutant_silent_struct_bound, mutant_struct_file, false /* score only */ );
	mutant_silent_data.write_silent_struct( mutant_silent_struct_unbound, mutant_struct_file, false /* score only */ );

	///// TODO: JUST FOR TESTING ///////
	pose.dump_pdb( mutation_tag + "_bound" + ".pdb" );
	////////////////////////////////////

	// Now reset the fold tree of the unbound pose back to its original fold tree
	unbound_pose.fold_tree( saved_foldtree );

	// Reset the atom_pair_constraint weight to intitial
	sfxn->set_weight( atom_pair_constraint, initial_cst_wt );

	return complex_score;
}

/////////////////////////////////////////////////////////////////////////////////
// This just makes sure csts are off for scoring
core::Real get_complex_score(
	core::scoring::ScoreFunctionOP & sfxn,
	core::pose::Pose & pose )
{
	using namespace basic::options;
	using namespace core::scoring;

	core::Real initial_cst_wt = sfxn->get_weight( atom_pair_constraint );
	sfxn->set_weight( atom_pair_constraint, 0.0 );
	core::Real complex_score( 0.0 );
	complex_score = (*sfxn)( pose );
	sfxn->set_weight( atom_pair_constraint, initial_cst_wt );
	// TODO: later might want to actually write pose to silent file here

	return complex_score;
}

/////////////////////////////////////////////////////////////////////////////////
core::Real calculate_binding_energy(
	core::scoring::ScoreFunctionOP & sfxn,
	core::pose::Pose & pose )
{
	utility::vector1< std::pair< core::Size, char > > empty_mutations;
	return calculate_binding_energy( sfxn, pose, empty_mutations );
}

/////////////////////////////////////////////////////////////////////////////////
core::Real get_unbound_protein_score(
	core::scoring::ScoreFunctionOP & sfxn,
	core::pose::Pose const & pose )
{

	using namespace basic::options;
	using namespace core::scoring;

	// TODO: Want to do packing and minimization just of the interface residues?
	// For now just pack and minimize the whole structure
	// If want to calculate more accurate unbound protein score, then repeat this many times
	// and take the average score
	// If there aren't any protein mutations and don't actually care about the unbound protein
	// score then protein_pack_reps can be set to 0
	core::Size pack_reps( option[ protein_pack_reps ]() );
	if ( pack_reps < 1 ) {
		return 0.0;
	}
	utility::vector1< core::Real > unbound_protein_scores;
	core::Real total_unbound_protein( 0.0 );
	// Get the starting protein pose
	core::pose::Pose protein_pose = pose;

	core::kinematics::FoldTree separation_foldtree;

	separation_foldtree = get_simple_RNA_protein_fold_tree( protein_pose );

	// Set the fold tree of the protein pose to the new separation fold tree
	// Should prevent problems after deletion of RNA residues
	protein_pose.fold_tree( separation_foldtree );

	delete_RNA_from_pose( protein_pose );

	for ( core::Size i=1; i<= pack_reps; ++i ) { // Default 10 reps

		core::pose::Pose protein_pose_copy = protein_pose;
		unbound_protein_minimization( protein_pose_copy, sfxn );
		if ( !option[min_only]() ) {
			simple_packing( protein_pose_copy, sfxn );
			unbound_protein_minimization( protein_pose_copy, sfxn );
		}

		// Scores of individual components
		core::Real unbound_protein_score;

		unbound_protein_score = (*sfxn)( protein_pose_copy );
		TR << "Unbound protein score: " << unbound_protein_score << std::endl;
		// append the protein score to the list
		unbound_protein_scores.push_back( unbound_protein_score );
		total_unbound_protein+=unbound_protein_score;
	}
	// Find the average unbound protein score
	core::Real avg_unbound_protein_score( total_unbound_protein / pack_reps );
	// TODO: Find and print the standard deviation (?)

	return avg_unbound_protein_score;

}

/////////////////////////////////////////////////////////////////////////////////
void mutate_and_score_RNP() {
	using namespace core::scoring;
	using namespace basic::options;
	using namespace core::import_pose::pose_stream;

	// Get the wildtype structure
	core::pose::Pose wt_pose;

	core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// This function is supposed to be smart about setting up the fold tree, I am using it
	// especially so that the jumps will get set up nicely so we don't see weird things in
	// minimization
	wt_pose = *(core::import_pose::initialize_pose_and_other_poses_from_command_line( rsd_set ));

	wt_pose.fold_tree().show( TR.Debug );

	// Check that the fold tree has at least one jump in it, if not, apply fold tree that
	// separates the RNA and protein by a jump (so torsion changes in one don't propagate
	// to the other) and also separates RNA and protein residues into different chains
	// based on ordering in the PDB
	// ( i.e. example 1 PDB order: RNA-protein-RNA --> RNA chain A - protein chain B - RNA chain C i
	// and jump 1 separates all RNA and protein;
	// example 2 PDB order: RNA-protein --> RNA chain A - protein chain B )

	if ( wt_pose.fold_tree().num_jump() < 1 ) {
		TR << "Default fold tree does not contain any jumps, which probably means you input a PDB that does not contain separate chains." << std::endl;
		TR << " Will try to guess a fold tree that separates RNA and protein into different chains!" << std::endl;
		core::kinematics::FoldTree new_fold_tree;
		new_fold_tree = get_simple_RNA_protein_fold_tree( wt_pose );
		wt_pose.fold_tree( new_fold_tree );
	}


	// Set up constraints if a constraint file was provided
	// These will only be turned on during minimization for now
	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		core::scoring::constraints::ConstraintSetOP cstset( new core::scoring::constraints::ConstraintSet() );
		cstset = core::scoring::constraints::ConstraintIO::read_constraints( option[OptionKeys::constraints::cst_file][1], cstset, wt_pose);
		wt_pose.constraint_set(cstset);
	}

	// Set up the score function
	ScoreFunctionOP sfxn;
	sfxn = setup_score_function();

	sfxn->show( wt_pose );

	// Set up the mutations
	utility::vector1< utility::vector1< std::pair< core::Size, char > > > mutations;
	utility::vector1< std::string > input_mutant_structs;
	if ( basic::options::option[ mutfile ].user() ) {
		std::pair< utility::vector1< utility::vector1< std::pair< core::Size, char > > >, utility::vector1< std::string > > mut_input;
		mut_input = get_mutations_from_file( basic::options::option[ mutfile ](), wt_pose );
		mutations = mut_input.first;
		input_mutant_structs = mut_input.second;
	} else { // Just wildtype
		TR << "Just running the wildtype structure..." << std::endl;
		TR << "Doing extra minimization and packing for the wildtype structure" << std::endl;
		sfxn->show( wt_pose );
		simple_minimization( wt_pose, sfxn );
		sfxn->show( wt_pose );
		if ( !option[ min_only ]() ) {
			simple_packing( wt_pose, sfxn );
			sfxn->show( wt_pose );
			simple_minimization( wt_pose, sfxn );
			sfxn->show( wt_pose );
		}
		calculate_binding_energy( sfxn, wt_pose );
	}

	utility::vector1< std::pair< core::Real, core::pose::Pose > > all_scores_and_poses;

	// Loop through the mutations
	for ( core::Size mutant_index = 1; mutant_index <= mutations.size(); ++mutant_index ) {
		core::pose::Pose mutant_pose, starting_mutant_pose;
		if ( (input_mutant_structs.size() > 0) and input_mutant_structs[mutant_index] != "none" ) {
			utility::vector1< core::Real > complex_scores;
			core::Real min_complex_score( 1000000000.0 );
			// Figure out which residues were rebuilt externally (listed in REBUILD_RESIDUES.txt file)
			// Get the directory from the file name
			core::Size last = 0;
			core::Size next = 0;
			std::string delimiter = "/";
			while ( (next = input_mutant_structs[mutant_index].find(delimiter,last)) != std::string::npos ) {
				last = next + delimiter.length();
			}
			std::string path = input_mutant_structs[mutant_index].substr(0,last);

			utility::vector1< core::Size > residues;
			std::ifstream residue_input;
			residue_input.open( (path + "REBUILD_RESIDUES.txt").c_str() );
			core::Size residue;
			while ( residue_input >> residue ) {
				residue += basic::options::option[ prot_offset_num ]();
				residue += basic::options::option[ RNA_offset_num ]();
				// The RNA residue numbering is 1 based in REBUILD_RESIDUES.txt
				// vs 0 based in the mutfile, so RNA_offset_num is off by 1
				// TODO: This will almost certainly cause problems
				residue -= 1;
				residues.push_back(residue);
			}

			TR << "These residues were rebuilt externally, and they and their neighboring residues will be rebuilt:" << std::endl;
			TR << residues << std::endl;
			TR << "These residues correspond to: " << std::endl;
			for ( core::Size i=1; i<=residues.size(); ++i ) {
				TR << "WT residue " << residues[i] << wt_pose.residue( residues[i] ).name1() << std::endl;
			}

			// Also add any residues that we're going to mutate
			for ( core::Size r=1; r<=mutations[mutant_index].size(); ++r ) {
				residues.push_back( mutations[ mutant_index ][r].first );
			}

			utility::vector1< core::Size > neighbor_residues =
				find_residues_around_mutation( wt_pose, residues );

			// Then we should read in the silent file and loop through its structures
			PoseInputStreamOP input;
			input = PoseInputStreamOP( new SilentFilePoseInputStream( input_mutant_structs[mutant_index] ));
			while ( input->has_another_pose() ) {
				input->fill_pose( starting_mutant_pose, *rsd_set );
				// Loop through mutation protocol many times
				core::Size nreps = option[ Nreps ](); // Default=10
				for ( core::Size n=1; n<=nreps; ++n ) {
					mutant_pose = starting_mutant_pose;
					// Make the mutations, if there are any
					if ( mutations[ mutant_index ].size() > 0 || neighbor_residues.size() > 0 ) {
						make_mutations( mutant_pose, sfxn, mutations[ mutant_index ], neighbor_residues, true );
					}
					core::Real complex_score = get_complex_score( sfxn, mutant_pose );
					TR << "Score of the complex: " << complex_scores << std::endl;
					complex_scores.push_back( complex_score );
					if ( complex_score < min_complex_score ) {
						min_complex_score = complex_score;
					}
					all_scores_and_poses.push_back( std::make_pair( complex_score, mutant_pose ) );
				}
			}
			TR << "MINIMUM COMPLEX SCORE: " << min_complex_score << std::endl;
			core::Real sum_complex_scores(0.0);
			for ( core::Size i = 1; i<=complex_scores.size(); ++i ) {
				sum_complex_scores+=complex_scores[i];
			}
			TR << "AVERAGE COMPLEX SCORE: " << sum_complex_scores/complex_scores.size() << std::endl;
		} else { // there is no input mutant stucture, so create it here
			utility::vector1< core::Size > neighbor_residues =
				find_residues_around_mutation( wt_pose, mutations[ mutant_index ] );
			// Generate a mutant structure and calculate the complex scores N times
			core::Size nreps = option[ Nreps ](); // Default=10
			for ( core::Size n=1; n<=nreps; ++n ) { //TODO: change to command line option?
				// Make a copy of the wildtype pose
				mutant_pose = wt_pose;
				// Make the mutations
				make_mutations( mutant_pose, sfxn, mutations[ mutant_index ], neighbor_residues, false );
				// Score the mutant structure and calculate binding energy
				core::Real complex_score = get_complex_score( sfxn, mutant_pose );
				TR << "Score of the complex: " << complex_score << std::endl;
				// Add the score and the structure to a vector
				all_scores_and_poses.push_back( std::make_pair( complex_score, mutant_pose ) );
			}
		}
		// Figure out the top 20% of the structures and scores
		// Output the average of the top 20% of scores
		utility::vector1< core::Real > all_scores;
		core::Real min_score = 1000000000.0;
		core::pose::Pose lowest_scoring_pose;
		for ( core::Size i =1; i<=all_scores_and_poses.size(); ++i ) {
			all_scores.push_back( all_scores_and_poses[i].first );
			if ( all_scores_and_poses[i].first < min_score ) {
				min_score = all_scores_and_poses[i].first;
				lowest_scoring_pose = all_scores_and_poses[i].second;
			}
		}
		std::sort( all_scores.begin(), all_scores.end() );
		core::Real sum_top_scores( 0.0 );
		core::Size num_poses_to_avg = 1;
		if ( core::Size(0.2*all_scores.size()) > num_poses_to_avg ) {
			num_poses_to_avg = core::Size(0.2*all_scores.size());
		}
		for ( core::Size i=1; i<=num_poses_to_avg; ++i ) {
			sum_top_scores += all_scores[i];
		}
		core::Real avg_top_scores = sum_top_scores/num_poses_to_avg;
		TR << "Average of top " << num_poses_to_avg << " poses: " << avg_top_scores << std::endl;
		// Output the lowest scoring structure
		std::string mut_tag = make_mutation_tag( mutations[ mutant_index ] );
		lowest_scoring_pose.dump_pdb("complex_lowest_scoring_" + mut_tag + ".pdb");

		// Figure out the unbound protein score
		core::Real avg_unbound_protein_score(0.0);
		core::Size pack_reps = option[ protein_pack_reps ]();
		if ( pack_reps > 0 ) {
			// For now use the lowest_scoring_pose to start unbound protein score calculation
			// TODO: might want to change this later
			avg_unbound_protein_score = get_unbound_protein_score( sfxn, lowest_scoring_pose );
			TR << "Average unbound protein score: " << avg_unbound_protein_score << std::endl;
		}
	}

}

int main( int argc, char ** argv ) {

	try {
		using namespace basic::options;
		NEW_OPT( mutfile, "File listing mutations", "test_mutfile" );
		NEW_OPT( move_backbone, "Move the backbone?", false );
		NEW_OPT( move_protein_backbone, "Move the protein backbone?", false );
		NEW_OPT( min_jumps, "Minimize the jumps?", false ); // Added 1-12-16: before this jumps were minimized
		NEW_OPT( out_prefix, "Prefix for the output files", "" );
		NEW_OPT( unbound_protein, "Structure to use for the unbound protein", "" );
		NEW_OPT( unbound_RNA, "Structure to use for the unbound RNA", "" );
		NEW_OPT( unbound_wtRNA_ens, "An ensemble of unbound RNA structures for the wt sequence", "" );
		NEW_OPT( mutate_unbound_RNA, "Should the unbound RNA get mutated to the new sequence?", false );
		NEW_OPT( prot_offset_num, "Offset to apply for mutfile numbering, how many protein residues come before RNA to mutate?", 0 );
		NEW_OPT( RNA_offset_num, "Offset to apply for mutfile numbering for RNA, does RNA numbering start at 1?", 0 );
		NEW_OPT( protein_start_res, "protein start residue", 0 );
		NEW_OPT( protein_end_res, "protein end residue", 0 );
		NEW_OPT( relax_cutoff_dist, "Residues within this cutoff distance of mutations will be relaxed", 20);
		NEW_OPT( min_only, "Minimize only? (Don't re-pack)", false );
		NEW_OPT( dump_bound_rna, "Dump a pdb containing only the bound RNA (not minimized alone)", false);
		NEW_OPT( bound_rna_dump_tag, "Tag for the structure from dump_bound_rna", "");
		NEW_OPT( protein_pack_reps, "Number of times to repack the unbound protein structure then avg", 1);
		NEW_OPT( Nreps, "Number of times to calculate complex score (top 20% will be averaged)", 10);

		devel::init( argc, argv );
		mutate_and_score_RNP();
	} catch ( utility::excn::Exception const & e ) {
		TR << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
