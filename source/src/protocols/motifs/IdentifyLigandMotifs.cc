// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IdentifyLigandMotifs.cc
/// @brief .cc fils for IdentifyLigandMotifs protocol. Protocol object reads in a pdb file(s) and outputs motifs for protein-ligand interactions in .pdb and .motifs format. App originally written by mdsmith, optimized and converted to protocol by Ari Ginsparg.

// libRosetta headers
#include <basic/options/option.hh>

#include <sstream>
#include <string>
#include <core/types.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>

////////////////////////////
////////////////////////////
////////////////////////////
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoreType.hh>

///////////////////////////////
//new, from upstream filter
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.fwd.hh>
/////////////////////////////////////////

#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>
//#include <core/pack/types.hh>

#include <utility/graph/Graph.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
//#include <core/id/AtomID_Map.Pose.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/scoring/mm/MMTorsionLibrary.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/util.hh>

#include <basic/basic.hh>

#include <basic/database/open.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

#include <protocols/dna/util.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/Motif.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>

#include <protocols/motifs/IdentifyLigandMotifs.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>

#include <protocols/motifs/motif_utils.hh>

// Time profiling header
#include <time.h>

using namespace core;
using namespace basic;
using namespace ObjexxFCL;
using namespace pose;
using namespace chemical;
using namespace scoring;
using namespace optimization;
using namespace options;
using namespace OptionKeys;

using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "apps.pilot.ligand_motifs", basic::t_info );

//default constructor
IdentifyLigandMotifs::IdentifyLigandMotifs()
{
	//set up values from options
	motif_pdb_output_path_ = option[ OptionKeys::motifs::ligand_motif_output_directory_name ];
	motif_file_output_ = option[ OptionKeys::motifs::ligand_motif_output_file_name ];
	//default true
	bool output_motifs_ = option[ OptionKeys::motifs::output_motifs ];
	//default true
	bool output_motifs_as_pdb_ = option[ OptionKeys::motifs::output_motifs_as_pdb ];

	TR << "motif_pdb_output_path_: " << motif_pdb_output_path_ << std::endl;
	TR << "motif_file_output_: " << motif_file_output_ << std::endl;
	TR << "output_motifs_: " << output_motifs_ << std::endl;
	TR << "output_motifs_as_pdb_: " << output_motifs_as_pdb_ << std::endl;
}

//destructor
IdentifyLigandMotifs::~IdentifyLigandMotifs() = default;

//get hbond score of ligand-residue interaction
core::Real
IdentifyLigandMotifs::get_hbond_score(
	core::pose::Pose & pose,
	core::Size pos1,
	core::Size pos2,
	core::scoring::ScoreFunction & scorefxn
)
{
	core::scoring::EnergyMap hbond_map;
	//run score function on pose
	//moving this so scoring doesn't occur every time this is called (likely only has to be done once)
	//scorefxn(pose);
	//context independent first, then dependent
	scorefxn.eval_ci_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, hbond_map );
	scorefxn.eval_cd_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, hbond_map );
	return hbond_map[ core::scoring::hbond_sc ];
}

//write a good motif interaction to a MotifLibrary
void
IdentifyLigandMotifs::output_single_motif_to_MotifLibrary(
	core::pose::Pose & src_pose,
	protocols::motifs::MotifLibrary & motifs,
	std::string & pdb_name,
	core::Size prot_pos,
	core::Size lig_pos,
	utility::vector1< core::Size > & lig_atoms,
	core::Real pack_score,
	core::Real hb_score
)
{
	core::Size prot_pos_pdb( src_pose.pdb_info()->number( prot_pos ) );

	conformation::Residue const & protres( src_pose.residue( prot_pos ) );
	conformation::Residue const & ligres( src_pose.residue( lig_pos ) );

	//Here is where we deal with the 3 ligand motif atoms
	//shouldn't have to do this in making the motif for the library?
	/*
	for ( core::Size atom_i = 1; atom_i <= ligres.natoms(); ++atom_i ) {

	if ( atom_i != lig_atoms[1] && atom_i != lig_atoms[2] && atom_i != lig_atoms[3] ) {
	numeric::xyzVector< core::Real > newpos(ligres.xyz(atom_i));
	newpos[1] = newpos[1]+100;
	ligres.set_xyz(atom_i, newpos);
	}
	}
	*/

	protocols::motifs::Motif motif(protres, ligres, lig_atoms);

	//add a remark to the motif so that more is known about it and can  be called later
	//want to hold source pdb name, residue name, ligand name, and packing/hbond scores
	std::string pdb_name_str = pdb_name;

	//need to break up the pdb_name_str  so that it does not have any file path leading up to it; only take characters beyond the final /

	std::string build_string = "";

	//iterate through file name character by character
	for ( core::Size string_iterator = 0; string_iterator < pdb_name_str.size(); string_iterator++ ) {
		//append character to build string
		build_string += pdb_name_str[string_iterator];

		//reset build string if encounter a /
		if ( pdb_name_str[string_iterator] == '/' ) {
			build_string = "";
		}
	}

	pdb_name_str = build_string;

	std::string res_str = src_pose.residue_type( prot_pos ).name3() + string_of( prot_pos_pdb );
	std::string lig_str = src_pose.residue_type( lig_pos ).name3();
	std::string score_str = "Packing_score:" + std::to_string(pack_score) + "_Hbond_score:" + std::to_string(hb_score);

	motif.store_remark(pdb_name_str + "_" + res_str + "_" + lig_str + "_" + score_str);

	motifs.add_to_library( motif ); //Add motif to library
}

//write a good motif interaction out as a pdb file
void
IdentifyLigandMotifs::output_single_motif_to_pdb(
	core::pose::Pose & src_pose,
	std::string & pdb_name,
	core::Size prot_pos,
	core::Size lig_pos,
	utility::vector1< core::Size > & lig_atoms
)
{

	core::Size prot_pos_pdb( src_pose.pdb_info()->number( prot_pos ) );
	char prot_pos_chain( src_pose.pdb_info()->chain( prot_pos ) );

	//need to figure out how to change the output path to something customizable in the flags
	//std::string output_path( "Ligand_motif_dir/" );
	//output path needs to already exist or script will break here!
	//default will output to local directory called "Ligand_Motif_dir"
	std::string output_path(motif_pdb_output_path_);
	std::string delimiter( "_" );
	std::string extension( ".pdb" );
	std::string motif_file_name(
		output_path +
		src_pose.residue_type( prot_pos ).name3() + string_of( prot_pos_pdb ) + string_of( prot_pos_chain ) + delimiter );


	conformation::Residue protres( src_pose.residue( prot_pos ) );
	conformation::Residue ligres( src_pose.residue( lig_pos ) );
	//Here is where we deal with the 3 ligand motif atoms
	for ( core::Size atom_i = 1; atom_i <= ligres.natoms(); ++atom_i ) {

		if ( atom_i != lig_atoms[1] && atom_i != lig_atoms[2] && atom_i != lig_atoms[3] ) {
			numeric::xyzVector< core::Real > newpos(ligres.xyz(atom_i));
			newpos[1] = newpos[1]+100;
			ligres.set_xyz(atom_i, newpos);
		}
	}

	Pose pose;
	pose.append_residue_by_jump( protres, 1 );
	pose.append_residue_by_jump( ligres, 1 );
	pose.conformation().insert_chain_ending( 1 );

	core::Size lig_pos_pdb( src_pose.pdb_info()->number( lig_pos ) );
	char lig_pos_chain( src_pose.pdb_info()->chain( lig_pos ) );
	std::string lig_talk = "Ligatoms" + delimiter + string_of( lig_atoms[1] )  + delimiter + string_of( lig_atoms[2] )  + delimiter + string_of( lig_atoms[3]) ;
	std::string motif_file_name_2 = motif_file_name + src_pose.residue_type( lig_pos ).name1() + string_of( lig_pos_pdb ) +  delimiter + lig_talk + delimiter +  string_of( lig_pos_chain ) + delimiter;

	std::string motif_file_name_3 = motif_file_name_2 + pdb_name;
	//only add extension if not already ending in ".pdb" + extension;
	if ( motif_file_name_3[motif_file_name_3.length() - 4] != '.' && motif_file_name_3[motif_file_name_3.length() - 3] != 'p' && motif_file_name_3[motif_file_name_3.length() - 2] != 'd' && motif_file_name_3[motif_file_name_3.length() - 1] != 'b' ) {
		motif_file_name_3 = motif_file_name_3 + extension;
	}

	TR << "Writing " << motif_file_name_3 << std::endl;

	core::io::pdb::dump_pdb( pose, motif_file_name_3 );

}


//get packing score of ligand-residue interaction
core::Real IdentifyLigandMotifs::get_packing_score(
	core::pose::Pose & pose,
	core::Size pos1,
	core::Size pos2,
	core::scoring::ScoreFunction & scorefxn
)
{
	core::scoring::EnergyMap pack_map;
	scorefxn.eval_ci_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, pack_map );
	return pack_map[ core::scoring::fa_atr ];
}

//get atom trios
void IdentifyLigandMotifs::get_atom_trios(
	utility::vector1< utility::vector1< Size > > & motif_indices_list,
	utility::vector1< utility::vector1< Size > > & all_motif_indices,
	core::conformation::Residue const & ligres
)
{
	for ( core::Size atom_i = 1; atom_i <= ligres.natoms(); ++atom_i ) {
		if ( ligres.atom_type(atom_i).is_hydrogen() ) {
			break;
		}

		TR << "in atom iterate block, atom num is " << atom_i <<   std::endl;
		std::string const atom_name = ligres.atom_name(atom_i);
		TR << "atom name is " << atom_name <<  std::endl;

		// This is a for loop to iterate over each atom's connected atoms:
		//core::conformation::Residue::AtomIndices atom_i_connects(  ligres.bonded_neighbor( atom_i ) );
		//for ( core::Size atom_j = 1; atom_j <= atom_i_connects.size(); ++ atom_j ) {
		for ( core::Size atom_j: ligres.bonded_neighbor( atom_i ) ) {
			//if ( ligres.atom_type(atom_i_connects[atom_j]).is_hydrogen() ) {
			if ( ligres.atom_type(atom_j).is_hydrogen() ) {
				break;
			}
			TR << "ATOM j: " << atom_j << " Name: " << ligres.atom_name(atom_j) << std::endl;

			// This is the next for loop to find connects for the second atom, giving us the final atom number (atom k)
			//core::conformation::Residue::AtomIndices atom_j_connects(  ligres.bonded_neighbor( atom_i_connects[atom_j] ) );
			//for ( core::Size atom_k = 1; atom_k <= atom_j_connects.size(); ++ atom_k ) {
			for ( core::Size atom_k: ligres.bonded_neighbor( atom_j ) ) {
				if ( ligres.atom_type(atom_k).is_hydrogen() ) {
					break;
				}
				TR << "ATOM k: " << atom_k << " Name: " << ligres.atom_name(atom_k) << std::endl;

				chemical::AtomType const & atom_i_type(ligres.atom_type(atom_i));
				std::string atom_i_name = atom_i_type.atom_type_name();
				chemical::AtomType const & atom_j_type(ligres.atom_type(atom_j));
				std::string atom_j_name = atom_j_type.atom_type_name();
				chemical::AtomType const & atom_k_type(ligres.atom_type(atom_k));
				std::string atom_k_name = atom_k_type.atom_type_name();

				TR << "Connected triplet is: " << atom_i << ", type is " << atom_i_name  << "; ";
				TR << atom_j << ", type is " << atom_j_name << "; " ;
				TR << atom_k << ", type is " << atom_k_name << " " << std::endl;
				if ( atom_i != atom_k ) {

					//make the 3 atom vector
					utility::vector1< Size > cur_motif_indices;
					cur_motif_indices.push_back( atom_i );
					cur_motif_indices.push_back( atom_j );
					cur_motif_indices.push_back( atom_k );

					TR << "atom_i: " << atom_i << std::endl;
					TR << "atom_i to atom_j: " << atom_j << std::endl;
					TR << "atom_j to atom_k: " << atom_k << std::endl;

					//check if current index is a duplicate

					bool current_is_duplicate ( false );
					if ( !motif_indices_list.empty() ) {
						for ( core::Size cur_motif_check = 1; cur_motif_check <= motif_indices_list.size(); ++ cur_motif_check ) {
							// run a test to see if vectors are identical
							utility::vector1< Size > cur_mainlist ( motif_indices_list[cur_motif_check] );
							utility::vector1< Size > sort_from_mainlist (cur_mainlist) ;
							std::sort(  sort_from_mainlist.begin(),  sort_from_mainlist.end() );
							utility::vector1< Size > sort_from_curindex (cur_motif_indices) ;
							std::sort( sort_from_curindex.begin(), sort_from_curindex.end() );
							if ( sort_from_mainlist[1] == sort_from_curindex[1] && sort_from_mainlist[2] == sort_from_curindex[2] && sort_from_mainlist[3] == sort_from_curindex[3] ) {
								current_is_duplicate = true;
							}
						}
					}
					//Making master list to see what's getting pruned with my check
					all_motif_indices.push_back( cur_motif_indices );
					//if current index isn't a duplicate, add it to the vector
					if ( !current_is_duplicate ) {
						motif_indices_list.push_back( cur_motif_indices );
					}
				}
			}
		}
	}
}

//main function, iterate over ligand (broken into adjacent 3 atom trios) and identify how they interact with nearby residues
void
IdentifyLigandMotifs::process_for_motifs(
	Pose & pose,
	std::string & pdb_name,
	protocols::motifs::MotifLibrary & motifs
)
{

	core::scoring::ScoreFunctionOP scorefxn(ScoreFunctionFactory::create_score_function( "ligand.wts" ));

	//score the pose (in theory only has to be done once)
	scorefxn->score(pose);

	core::Size nres( pose.size() );

	// Loop over positions, skipping non-amino acid

	//////////////////////
	//////////////////////
	//////////////////////

	//ligand residue loop
	for ( core::Size lig_pos = 1 ; lig_pos <= nres ; ++lig_pos ) {
		ResidueType const & lig_type( pose.residue_type( lig_pos ) );

		// .is_ligand() comes from  src/core/chemical/ResidueType.cc (I think)
		//look at ligands only
		if (  !lig_type.is_ligand() ) continue;
		TR << "in ligand splitter block, found my ligand, lig_pos is " << lig_pos << std::endl;

		// This is to make a ligres object once we find our ligand
		conformation::Residue const & ligres( pose.residue( lig_pos ) );

		//we want to rule out collecting motifs from "small" ligands and ions
		//going to define "small" as fewer than 3 non-hydrogen and non-virtual atoms (hydrogen and virtual will not count at all)
		//if we do not have 3 heavy atoms, we can't make a proper motif
		core::Size num_virt_atoms = ligres.n_virtual_atoms();

		core::Size num_non_virt_heavy_atoms = ligres.nheavyatoms() - num_virt_atoms;

		//continue if we have 2 or fewer non-virtual atoms
		if ( num_non_virt_heavy_atoms < 3 ) {
			TR << "Ignoring ligand " << ligres.name() << " that has " << num_non_virt_heavy_atoms << " heavy and non-virtual atoms. It has " << num_virt_atoms << " virtual atoms and " << ligres.nheavyatoms() << " heavy atoms." << std::endl;
			continue;
		}

		//get ligand atom trio block
		//all unique trios
		utility::vector1< utility::vector1< Size > > motif_indices_list;
		//all adjacent trios (includes duplicates)
		utility::vector1< utility::vector1< Size > > all_motif_indices;

		get_atom_trios(motif_indices_list, all_motif_indices, ligres);

		TR << "Total 3 atoms in pruned indices list is: " << motif_indices_list.size()  << std::endl;
		TR << "Total 3 atoms in unpruned indices list is: " << all_motif_indices.size()  << std::endl;

		////////////////////////////
		////////////////////////////
		////////////////////////////
		//loop over amino acids to attempt to pair with ligand for identifying motifs
		for ( core::Size prot_pos = 1 ; prot_pos <= nres ; ++prot_pos ) {

			//break motif identification into its own function for better readability
			ligand_to_residue_analysis(lig_pos, prot_pos, pose, pdb_name, motifs, scorefxn, motif_indices_list);

		}

		//Here we're going to check to see what's in motif_indices_list
		for ( core::Size motif_position = 1; motif_position <= motif_indices_list.size(); ++ motif_position ) {
			utility::vector1< Size > cur_trip ( motif_indices_list[motif_position] );
			TR << "Motif index contains: " << cur_trip[1] << "-" << cur_trip[2] << "-" << cur_trip[3] << std::endl;
		}
	}
}

void
IdentifyLigandMotifs::ligand_to_residue_analysis(
	core::Size lig_pos,
	core::Size prot_pos,
	core::pose::Pose & pose,
	std::string & pdb_name,
	protocols::motifs::MotifLibrary & motifs,
	core::scoring::ScoreFunctionOP scorefxn,
	utility::vector1< utility::vector1< Size > > & motif_indices_list
)
{
	ResidueType const & prot_type( pose.residue_type( prot_pos ) );
	if (  !prot_type.is_protein() ) return;
	if ( pose.residue( prot_pos ).name3() == "GLY" ) return;

	// map will automatically sort the "contacts" with the lowest total_score at the front of map
	std::map< Real, Size > contacts;
	std::map< Real, Size > distance_sorter;

	ResidueType const & lig_type( pose.residue_type( lig_pos ) );

	if (  !lig_type.is_ligand() ) {
		return;
	}

	//skip if distance between residue and ligand nbr atoms is greater than 1.5 times the sum of their nbr radii

	core::Real lig_nbr_radius = pose.residue( lig_pos ).nbr_radius();
	core::Real res_nbr_radius = pose.residue( prot_pos ).nbr_radius();
	core::Real test_val = 1.5 * (lig_nbr_radius + res_nbr_radius);
	core::Real lig_res_nbr_distance = pose.residue( prot_pos ).nbr_atom_xyz().distance(pose.residue( lig_pos ).nbr_atom_xyz());
	if ( lig_res_nbr_distance > test_val ) {
		return;
	}


	//TR << "past lig_type statement" << std::endl;
	Real pack_score = get_packing_score( pose, lig_pos, prot_pos, *scorefxn );
	Real hb_score = get_hbond_score( pose, lig_pos, prot_pos, *scorefxn );
	//TR << "before get_packing_score call" << std::endl;

	//TR << "after get_packing_score call" << std::endl;
	Real total_score = pack_score + hb_score;
	//adding in total_score restrain to be at least -2.0
	if ( ( pack_score < -0.5 || hb_score < -0.3) && total_score < -1.0 ) { //Enter into motif checking loop--For each interacting residue:

		//TR << "In pack_score statement" << std::endl;

		contacts[total_score] = lig_pos;
		TR << "Residue " << prot_pos << " passed energy cut with pack score: " << pack_score << ", hbond score: " << hb_score << ", for a total score of: " << total_score << std::endl;

		for ( core::Size motif_position = 1; motif_position <= motif_indices_list.size(); ++ motif_position ) { //for each 3 atom triplet from ligand
			bool resi_trip_match( false );

			utility::vector1< Size > cur_trip ( motif_indices_list[motif_position] );

			Real closest_distance(5.0);

			for ( core::Size cur_trip_pos = 1; cur_trip_pos <= cur_trip.size(); ++ cur_trip_pos ) {   //for each atom in ligand triplet
				for ( core::Size residue_atom_number = 1; residue_atom_number <= pose.residue( prot_pos ).nheavyatoms(); ++ residue_atom_number ) {   //for each atom in current residue


					if ( cur_trip[cur_trip_pos] <= pose.residue( lig_pos ).natoms() ) {
						//This is the potentially broken line
						Real atom_atom_distance( pose.residue( lig_pos ).xyz( cur_trip[cur_trip_pos]  ).distance( pose.residue( prot_pos ).xyz( residue_atom_number ) ) );


						if ( atom_atom_distance < 4.0 ) {

							resi_trip_match = true;
							if ( atom_atom_distance < closest_distance ) {
								closest_distance=atom_atom_distance;
							}
						}
					} else {
						TR << "illegal value of cur_trip[cur_trip_pos]" << std::endl;
					}
				}
			}

			//skip motif if closest distance is beyond 5.0 (won't be good for much)
			if ( closest_distance == 5.0 ) {
				TR << "Skipping motif whose closest atom-atom distance is no closer than 5.0 angstroms" << std::endl;
				continue;
			}


			//here is the end of what we do with each triplet, this is where we add triplet-residue to vector of pairs
			if ( resi_trip_match ) {
				distance_sorter[closest_distance] = motif_position;
			}
		}
		//Here is where we find the closest few triplets for the current residue (we are in a protein residue loop)

		utility::vector1< Size > top_triplets;
		Size sorter_size( distance_sorter.size() );

		/////////////////////
		/////////////////////
		if ( sorter_size != 0 ) {
			Size min(2); //This is the max number of motifs to output
			// TR << "Minimum is " << min << std::endl;
			Size stop_after_min_count(1);
			////////////////////////////////////////////////////////////////////////////////////////////
			for ( std::map< Real, Size >::const_iterator it( distance_sorter.begin() ),  end( distance_sorter.end() ); it != end; ++it ) {
				if ( stop_after_min_count == min ) { break; }
				Size my_cur_trip( it-> second );
				top_triplets.push_back( my_cur_trip );
				++stop_after_min_count;
			}
			///////////////////////////////////////////////////////////////////////////////////////
		}
		///////////////////
		///////////////////

		TR << "Top triplets contains " << top_triplets.size() << " items." << std::endl << "Top triplets are: " ;
		for ( core::Size top_trip_pos = 1; top_trip_pos <= top_triplets.size(); ++ top_trip_pos ) {
			Size this_triplet_number(top_triplets[top_trip_pos]);
			utility::vector1< Size > this_triplet( motif_indices_list[this_triplet_number] );
			TR << "Size of top_triplets:  " << top_triplets.size() << std::endl;
			TR << "Size of this_triplet:  " << this_triplet.size() << std::endl;
			TR <<  this_triplet_number << ": " << this_triplet[1] << "-" << this_triplet[2] << "-" <<  this_triplet[3]  ;

			//output motifs to library and .motifs file
			if ( output_motifs_ ) {
				output_single_motif_to_MotifLibrary( pose, motifs, pdb_name, prot_pos, lig_pos, this_triplet, pack_score, hb_score );
			}
			//output motifs as pdb files
			if ( output_motifs_as_pdb_ ) {
				output_single_motif_to_pdb( pose, pdb_name, prot_pos, lig_pos, this_triplet );
			}

		}
	}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//process inputted pdb file(s)
void
IdentifyLigandMotifs::process_file_list()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	//using namespace core::io::pdb;


	using namespace basic::options;
	//using namespace basic::options::OptionKeys;
	using namespace optimization;

	ScoreFunction scorefxn;

	//core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );


	scorefxn.set_weight( fa_atr, 1.00 );
	scorefxn.set_weight( fa_rep, 1.00 );
	scorefxn.set_weight( fa_sol, 1.00 );
	scorefxn.set_weight( fa_dun, 1.00 );
	scorefxn.set_weight( fa_pair, 1.00 );
	scorefxn.set_weight( p_aa_pp, 1.00 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );


	//removing this since we may not necessarily use -s or -l for file inputs. Instead, going to use a custome string input (which also allows for increased customization)
	//utility::vector1< std::string > pdb_files( start_files() );

	//Here is where I create the MotifLibrary object which I populate with the motifs
	//will feed this created library into the global variable at the end of the protocol
	//protocols::motifs::MotifLibrary motifs;

	//new code for pose input (using Rocco's recommendation for using src/core/import_pose/pose_stream/MetaPoseInputStream.hh)
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

	//I don't think we need to specify a residue type set, but keeping this in case we do/want to
	//core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	while ( input.has_another_pose() ) {
		core::pose::Pose pose;
		//this would use a specified RTS
		//input.fill_pose( pose, *rsd_set );
		//this doesn't
		input.fill_pose( pose );
		//now we have our pose from file

		//derive a pdb name to use
		std::string pdb_name = input.get_last_pose_descriptor_string();
		//run process_for_motifs
		process_for_motifs( pose, pdb_name, motif_library_ );
	}

	//decoupling write of motifs to disk so that it may be called separately if desired
	/*
	protocols::motifs::MotifCOPs motifcops = motif_library_.library();
	//motifcops_ = motifs.library();

	//std::string filename( "AllMattMotifsFile.motifs" );
	//use option to use custom name for motif file output
	std::string filename( motif_file_output_ );

	//if there is not period in the file name, append ".motifs" to the end
	if ( filename.find('.') == std::string::npos ) {
	filename = filename + ".motifs";
	}

	utility::io::ozstream motif_output_file( filename );

	for ( protocols::motifs::MotifCOP & motifcop: motifcops ) {
	motif_output_file << *motifcop;
	}
	*/
}

// @brief returns the motif_library_ within the class, in case it is needed for additional usage beyond the scope of this protocol
protocols::motifs::MotifLibrary
IdentifyLigandMotifs::get_motif_library()
{
	return motif_library_;
}

// @brief writes the motif_library_ to a file, using the motif_file_output_ as the file name prefix
void
IdentifyLigandMotifs::write_motifs_to_disk()
{
	protocols::motifs::write_motifs_to_disk(motif_library_, motif_file_output_);
}
