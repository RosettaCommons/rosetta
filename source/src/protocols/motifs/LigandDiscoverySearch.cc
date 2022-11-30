// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief aginsparg, ipatel, sthyme


//#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/io/pdb/pdb_writer.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/dna/setup.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/motifs/LigandMotifSearch.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/motifs/Motif.hh>
#include <core/chemical/residue_io.hh>
#include <protocols/motifs/MotifHit.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/motifs/BuildPosition.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/pose/PDBInfo.hh>

#include <core/chemical/GlobalResidueTypeSet.hh>

//#include <protocols/motifs/MotifLigandPacker.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/util.hh>
//#include <protocols/dna/DnaInterfacePacker.hh>
//#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/motifs/motif_utils.hh>
using namespace protocols::dna;

#include <basic/prof.hh>
#include <basic/Tracer.hh>
static basic::Tracer TR( "protocols.ligand_discovery.LigandDiscoverySearch" );

#include <core/import_pose/import_pose.hh> // Need since refactor

#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
using utility::string_split;

// c++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <functional>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

//this may solve atom clash issue
#include <core/pose/xyzStripeHashPose.hh>
#include <numeric/geometry/hashing/xyzStripeHash.hh>
#include <numeric/geometry/hashing/xyzStripeHash.fwd.hh>

#include <numeric/xyzVector.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.hh>
#include <protocols/ligand_docking/ligand_scores.hh>

#include <protocols/qsar/scoring_grid/AtrGrid.hh>
#include <protocols/qsar/scoring_grid/RepGrid.hh>
#include <protocols/qsar/scoring_grid/VdwGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbdGrid.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>

#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/ligand_dock_impl.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

#include <ObjexxFCL/FArray1D.hh>

#include <protocols/motifs/LigandDiscoverySearch.hh>

// Time profiling header
#include <time.h>

//for data type debugging
#include <typeinfo>

using namespace core;
using namespace basic;
using namespace chemical;
using namespace pack;
using namespace task;
using namespace scoring;
using namespace options;
using namespace OptionKeys;

////////////////////////////////////////////////////////////////////////////

//merge sort code as well as a comparator struct to be used with a standard heap

//used for sorting poses by ddg score to obtain best
std::vector<std::tuple<core::Real, core::pose::Pose, std::string>> merge(std::vector<std::tuple<core::Real, core::pose::Pose, std::string>> left, std::vector<std::tuple<core::Real, core::pose::Pose, std::string>> right)
{
	std::vector<std::tuple<core::Real, core::pose::Pose, std::string>> result;
	while ( (int)left.size() > 0 || (int)right.size() > 0 ) {
		if ( (int)left.size() > 0 && (int)right.size() > 0 ) {
			if ( std::get<0>(left.front()) <= std::get<0>(right.front()) ) {
				result.push_back(left.front());
				left.erase(left.begin());
			} else {
				result.push_back(right.front());
				right.erase(right.begin());
			}
		}  else if ( (int)left.size() > 0 ) {
			for ( int i = 0; i < (int)left.size(); i++ ) {
				result.push_back(left[i]);
			}
			break;
		}  else if ( (int)right.size() > 0 ) {
			for ( int i = 0; i < (int)right.size(); i++ ) {
				result.push_back(right[i]);
			}
			break;
		}
	}
	return result;
}

std::vector<std::tuple<core::Real, core::pose::Pose, std::string>> mergeSort(std::vector<std::tuple<core::Real, core::pose::Pose, std::string>> m)
{
	if ( m.size() <= 1 ) {
		return m;
	}

	std::vector<std::tuple<core::Real, core::pose::Pose, std::string>> left, right, result;
	int middle = ((int)m.size()+ 1) / 2;

	for ( int i = 0; i < middle; i++ ) {
		left.push_back(m[i]);
	}

	for ( int i = middle; i < (int)m.size(); i++ ) {
		right.push_back(m[i]);
	}

	left = mergeSort(left);
	right = mergeSort(right);
	result = merge(left, right);

	return result;
}

struct comparator
{
	bool operator()(std::tuple<core::Real, core::pose::Pose, std::string> const& a, std::tuple<core::Real, core::pose::Pose, std::string> const& b) const
	{
		return std::get<0>(a) < std::get<0>(b);
	}
};

////////////////////////////////////////////////////////////////////////////
//functions/constructors to set up protocol
//default constructor
//will need to use class functions to seed values for input pose, motif library, and ligand library
LigandDiscoverySearch::LigandDiscoverySearch()
{

}

//destructor
LigandDiscoverySearch::~LigandDiscoverySearch() = default;

//parameterized constructor to load in motif library, pdb, and ligand library
LigandDiscoverySearch::LigandDiscoverySearch(core::pose::PoseOP pose_from_PDB, protocols::motifs::MotifCOPs motif_library, utility::vector1<core::conformation::ResidueOP> all_residues, core::Size working_position)
{

	working_pose_ = pose_from_PDB;
	motif_library_ = motif_library;
	all_residues_ = all_residues;
	working_position_ = working_position;

	//derive motif_library_for_select_residue_ from motif_library and residue in working_pose_ and index working_position_
	motif_library_for_select_residue_ = get_motif_sublibrary_by_aa(working_pose_->residue(working_position_).name3());

}

//function to load in a library for the search protocol
void LigandDiscoverySearch::set_motif_library(protocols::motifs::MotifCOPs motif_library)
{
	motif_library_ = motif_library;
}

//function to load in a library for the search protocol; read in MotifLibrary object and convert to MotifCOPs
void LigandDiscoverySearch::set_motif_library(protocols::motifs::MotifLibrary motif_library)
{
	motif_library_ = motif_library.library();
}

//function to define the residue index that we will use for applying motifs for ligand placements
void LigandDiscoverySearch::set_working_position(core::Size working_position)
{
	working_position_ = working_position;
}

//return contents of working_position_
core::Size LigandDiscoverySearch::get_working_position()
{
	return working_position_;
}

//return contents of motif_library_
protocols::motifs::MotifCOPs LigandDiscoverySearch::get_motif_library()
{
	return motif_library_;
}

//function to load in a pose for the receptor
void LigandDiscoverySearch::set_working_pose(core::pose::PoseOP pose_from_PDB)
{
	working_pose_ = pose_from_PDB;
}

//function to load in a pose for the receptor, use a regular pose and convert to pointer
void LigandDiscoverySearch::set_working_pose(core::pose::Pose pose_from_PDB)
{
	pose::Pose pose_helper = pose_from_PDB;
	pose::PoseOP pose_pointer_helper(new pose::Pose(pose_helper));
	working_pose_ = pose_pointer_helper;
}

//function to get working_pose_
core::pose::PoseOP LigandDiscoverySearch::get_working_pose()
{
	return working_pose_;
}

//function to set ligand residues in all_residues_
void LigandDiscoverySearch::set_all_residues(utility::vector1<core::conformation::ResidueOP> all_residues)
{
	all_residues_ = all_residues;
}

//function to get all ligand residues
utility::vector1<core::conformation::ResidueOP> LigandDiscoverySearch::get_all_residues()
{
	return all_residues_;
}

//main function to run ligand discovery operations
//needs to have values set for working_pose_, motif_library_, and all_residues_
//parameter is a string to be a prefix name to use for outputted file names
core::Size LigandDiscoverySearch::discover(std::string output_prefix)
{
	//create tracer to identify points of the run
	static basic::Tracer ms_tr( "LigandDiscoverySearch_out", basic::t_info );

	//determine whether to kill due to bad initialization
	bool kill_bad_init = false;

	//kill if working position is out of bounds
	if ( working_position_ < 1 || working_position_ > working_pose_->size() ) {
		ms_tr << "Working position of " << working_position_ << " is invalid as it is not a valid index to access a residue in your pose of size " << working_pose_->size() << std::endl;
		kill_bad_init = true;
	}


	std::string discovery_position_residue = working_pose_->residue(working_position_).name3();

	//get motif sublibrary
	motif_library_for_select_residue_ = get_motif_sublibrary_by_aa(discovery_position_residue);

	if ( motif_library_for_select_residue_.size() == 0 ) {
		ms_tr << "We have no motifs to work with here. Exiting function." << std::endl;
		kill_bad_init = true;
	}

	if ( kill_bad_init == true ) {
		ms_tr << "We have at least 1 bad initial input, killing the attempt now." << std::endl;
		return -1;
	}

	//seed score functions
	//seed_score_functions();

	core::scoring::ScoreFunctionOP fa_rep_fxn_(new core::scoring::ScoreFunction());
	core::scoring::ScoreFunctionOP fa_atr_fxn_(new core::scoring::ScoreFunction());
	core::scoring::ScoreFunctionOP fa_atr_rep_fxn_(new core::scoring::ScoreFunction());
	core::scoring::ScoreFunctionOP whole_score_fxn_(ScoreFunctionFactory::create_score_function( "ligand.wts" ));
	core::scoring::ScoreFunctionOP whole_score_fxn_no_constraints_(ScoreFunctionFactory::create_score_function( "ligand.wts" ));

	//fa_rep_fxn_ = new core::scoring::ScoreFunction();
	fa_rep_fxn_->set_weight(core::scoring::fa_rep, 0.44);

	//fa_atr_fxn_ = new core::scoring::ScoreFunction();
	fa_atr_fxn_->set_weight(core::scoring::fa_atr, 0.44);

	//fa_atr_rep_fxn_ = new core::scoring::ScoreFunction();
	fa_atr_rep_fxn_->set_weight(core::scoring::fa_atr, 0.44);
	fa_atr_rep_fxn_->set_weight(core::scoring::fa_rep, 0.44);

	//whole_score_fxn_ = ScoreFunctionFactory::create_score_function( "ligand.wts" );
	whole_score_fxn_->set_weight(core::scoring::fa_intra_rep, 0.004);
	whole_score_fxn_->set_weight(core::scoring::fa_elec, 0.42);
	whole_score_fxn_->set_weight(core::scoring::hbond_bb_sc, 1.3);
	whole_score_fxn_->set_weight(core::scoring::hbond_sc, 1.3);
	whole_score_fxn_->set_weight(core::scoring::rama, 0.2);
	whole_score_fxn_->set_weight(core::scoring::chainbreak, 1.0);
	whole_score_fxn_->set_weight(core::scoring::coordinate_constraint, 1.0);
	whole_score_fxn_->set_weight(core::scoring::atom_pair_constraint, 1.0);
	whole_score_fxn_->set_weight(core::scoring::angle_constraint, 1.0);
	whole_score_fxn_->set_weight(core::scoring::dihedral_constraint, 1.0);
	whole_score_fxn_->set_weight(core::scoring::res_type_constraint, 1.0);

	//whole_score_fxn_no_constraints_ = ScoreFunctionFactory::create_score_function( "ligand.wts" );
	whole_score_fxn_no_constraints_->set_weight(core::scoring::fa_intra_rep, 0.004);
	whole_score_fxn_no_constraints_->set_weight(core::scoring::fa_elec, 0.42);
	whole_score_fxn_no_constraints_->set_weight(core::scoring::hbond_bb_sc, 1.3);
	whole_score_fxn_no_constraints_->set_weight(core::scoring::hbond_sc, 1.3);
	whole_score_fxn_no_constraints_->set_weight(core::scoring::rama, 0.2);
	whole_score_fxn_no_constraints_->set_weight(core::scoring::coordinate_constraint, 0);

	core::Size x_shift = 0;
	core::Size y_shift = 0;
	core::Size z_shift = 0;
	int x_bound_int = 0;
	int y_bound_int = 0;
	int z_bound_int = 0;

	//create protein atom matrix
	create_protein_representation_matrix(x_shift, y_shift, z_shift, x_bound_int, y_bound_int, z_bound_int);

	//seeding of rep and atr values with default cuttoffs
	core::Real fa_rep_cutoff = 100;
	core::Real fa_atr_cutoff = -5;


	if ( option[ OptionKeys::motifs::fa_rep_cutoff ].user() ) {
		ms_tr << "Using user-inputted fa_rep cutoff of: "  << option[ OptionKeys::motifs::fa_rep_cutoff ] << std::endl;
		fa_rep_cutoff = option[ OptionKeys::motifs::fa_rep_cutoff ];
	}

	if ( option[ OptionKeys::motifs::fa_atr_cutoff ].user() ) {
		ms_tr << "Using user-inputted fa_atr cutoff of: "  << option[ OptionKeys::motifs::fa_atr_cutoff ] << std::endl;
		fa_atr_cutoff = option[ OptionKeys::motifs::fa_atr_cutoff ];
	}


	//create vector to hold the top X placements
	comparator comparator_v = comparator();
	//this gets used unless the value for the number of best placements to collect is 0 (in which all are collected)
	std::vector < std::tuple<core::Real, core::pose::Pose, std::string>> best_placements;
	//dynamic cutoff to determine what pplacements to consider; should get more strict as we get placements
	//get ddg cutoff
	core::Real ddg_cutoff_value = 100000;
	if ( option[ OptionKeys::motifs::ddg_cutoff ].user() ) {
		ddg_cutoff_value = option[ OptionKeys::motifs::ddg_cutoff ];
	}
	core::Real score_cutoff = ddg_cutoff_value;

	//determine how many output files to keep
	core::Size best_pdbs_to_keep = 0;
	if ( option[ OptionKeys::motifs::best_pdbs_to_keep ].user() ) {
		best_pdbs_to_keep = option[ OptionKeys::motifs::best_pdbs_to_keep ];
	}

	ms_tr << "Starting to iterate through all ligands" << std::endl;

	//hold the number of placements that pass all filters and could enter the top 100 placements
	int passed_placement_counter = 0;

	//now we have the filtered motif library to work with, run through each  atom trio in each ligand and try to match it against all motifs for the residue
	for ( core::Size tracker = 1; tracker <= all_residues_.size(); ++tracker ) {

		bool ligand_added = false;

		//core::chemical::ResidueTypeCOP ligres(ref);
		//convert ligres to be a ResidueOP type
		core::conformation::ResidueOP ligresOP = all_residues_[tracker];

		ms_tr << "On ligand " << ligresOP->name() << std::endl;

		const core::Real lig_nbr_radius = ligresOP->nbr_radius();

		ms_tr << "NBR_RADIUS of ligand is: " << lig_nbr_radius << std::endl;

		int ligand_passing_counter  = 0;
		int ligand_clashing_counter = 0;

		//derive a mutable residue type, used in high res docker
		core::chemical::MutableResidueTypeOP lig_mrt( new core::chemical::MutableResidueType( ligresOP->type() ) );

		// This is to make an atomtypeset to get atomtype integers
		core::chemical::AtomTypeSetCOP atset = core::chemical::ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );

		//3D vector to hold data of all connected atom trios
		//contents of vector are as follows:
		//ligand_atom_trios[trio identifier: 1-#trios][atom identifier within trio: 1-3][trio atom metadata]
		//metadata consists of index (position 1) and numerical atom_type_index (position 2)

		utility::vector1<utility::vector1< utility::vector1< core::Size > >> ligand_atom_trios;

		//find all atom trios (that do not contain hydrogen) in the ligand
		ms_tr << "Finding all atom trios" << std::endl;
		for ( core::Size atom_i = 1; atom_i <= ligresOP->natoms(); ++atom_i ) {

			//utility::vector1< core::Size > atom_vector;
			if ( ligresOP->atom_is_hydrogen(atom_i) ) { continue; }
			// This is a for loop to iterate over each atom's connected atoms:
			core::conformation::Residue::AtomIndices atom_i_connects(  ligresOP->bonded_neighbor( atom_i ) );
			for ( core::Size atom_j = 1; atom_j <= atom_i_connects.size(); ++ atom_j ) {
				if ( ligresOP->atom_is_hydrogen(atom_i_connects[atom_j]) ) { continue; }
				// This is the next for loop to find connects for the second atom, giving us the final atom number (atom k)
				core::conformation::Residue::AtomIndices atom_j_connects(  ligresOP->bonded_neighbor( atom_i_connects[atom_j] ) );
				for ( core::Size atom_k = 1; atom_k <= atom_j_connects.size(); ++ atom_k ) {
					if ( ligresOP->atom_is_hydrogen(atom_j_connects[atom_k]) ) { continue; }
					chemical::AtomType atom_i_type(ligresOP->atom_type(atom_i));
					if ( atom_i != atom_j_connects[atom_k] ) {

						//make the 3 atom vector
						utility::vector1< utility::vector1< core::Size > > cur_motif_indices;

						utility::vector1< core::Size > atom_i_vector;
						atom_i_vector.push_back( atom_i );
						atom_i_vector.push_back( atset->atom_type_index( ligresOP->atom_type(atom_i).atom_type_name() ) );

						utility::vector1< core::Size > atom_j_vector;
						atom_j_vector.push_back( atom_i_connects[atom_j] );
						atom_j_vector.push_back( atset->atom_type_index( ligresOP->atom_type(atom_i_connects[atom_j]).atom_type_name() ) );

						utility::vector1< core::Size > atom_k_vector;
						atom_k_vector.push_back( atom_j_connects[atom_k] );
						atom_k_vector.push_back( atset->atom_type_index( ligresOP->atom_type(atom_j_connects[atom_k]).atom_type_name() ) );

						cur_motif_indices.push_back( atom_i_vector);
						cur_motif_indices.push_back( atom_j_vector);
						cur_motif_indices.push_back( atom_k_vector);

						ligand_atom_trios.push_back(cur_motif_indices);
					}
				}
			}
		}

		//mini pose to represent a smaller region of the protein near the binding pocket for early energy calculations
		core::pose::PoseOP minipose(new pose::Pose);

		//now we have all trios, iterate through each trio and the motif sub-library
		ms_tr << "Looking through all atom trios" << std::endl;
		ms_tr << "#trios = " << ligand_atom_trios.size() << std::endl;

		//core::Size ligand_nbr_atom_index = ligresOP->nbr_atom();

		//iterate through each atom trio and try to use motifs to place the ligand
		for ( core::Size i = 1; i <= ligand_atom_trios.size(); ++i ) {

			ms_tr << "On trio # " << i << std::endl;
			ms_tr << "Trio is " << ligresOP->atom_name(ligand_atom_trios[i][1][1]) << " " << ligresOP->atom_name(ligand_atom_trios[i][2][1]) << " " << ligresOP->atom_name(ligand_atom_trios[i][3][1]) << std::endl;
			core::Size  trip_atom_1(ligand_atom_trios[i][1][1]);
			core::Size  trip_atom_2(ligand_atom_trios[i][2][1]);
			core::Size  trip_atom_3(ligand_atom_trios[i][3][1]);

			int clashing_counter = 0;
			int passing_counter = 0;
			int pose_atom_check_counter = 0;

			//run through the motif sublibrary and attempt to pair the ligand trio to the target residue based  on existing motifs
			//iterate through motif sublibrary
			ms_tr << "Looking through all motifs" << std::endl;
			ms_tr << "#motifs = " << motif_library_for_select_residue_.size() << std::endl;

			int motif_counter = 0;

			for ( auto motifcop : motif_library_for_select_residue_ ) {
				++motif_counter;

				if ( motif_counter % 10000 == 0 ) {
					ms_tr << "On motif #" << motif_counter <<std::endl;
				}
				//compare the atoms in  the ligand side  of the motif to the atom trio; continue if not a match
				//no need to check inverse order of trio, will get checked later

				//fixing error with comparison of signed and unsigned ints
				//making the motifcop res atoms be of size data type
				core::Size motif_ligand_atom_1 = motifcop->res2_atom1_int();
				core::Size motif_ligand_atom_2 = motifcop->res2_atom2_int();
				core::Size motif_ligand_atom_3 = motifcop->res2_atom3_int();

				//if motif atoms don't match trio atoms, continue
				if ( motif_ligand_atom_1 != ligand_atom_trios[i][1][2] || motif_ligand_atom_2 != ligand_atom_trios[i][2][2] || motif_ligand_atom_3 != ligand_atom_trios[i][3][2] ) {
					//at least one atom doesn't match, continue to next motif
					continue;
				}

				++pose_atom_check_counter;

				//use motif code to place the residue based on the pose residue
				motifcop->place_residue( working_pose_->residue(working_position_), *ligresOP, trip_atom_1, trip_atom_2, trip_atom_3 , true );

				//bool to determine if placed residue clashes against the backbone
				bool has_clashing = false;

				//check if the placement clashes
				has_clashing = ligand_clash_check(ligresOP, x_shift, y_shift, z_shift, x_bound_int, y_bound_int, z_bound_int);

				//debugging code for clash checking
				/*
				std::string pdb_name_early = output_prefix + "_ResPos_" + std::to_string(working_position_) + "_ResID_" + discovery_position_residue + "_Trio" + std::to_string(i) + "_" + ligresOP->name() + "_motif_" + motifcop->remark() + "_clashing_" + std::to_string(has_clashing) + ".pdb";
				working_pose_->append_residue_by_jump(*ligresOP, 1);
				core::io::pdb::dump_pdb(*working_pose_, pdb_name_early);
				working_pose_->delete_residue_slow(working_pose_->size());
				*/

				//continue because we clash
				if ( has_clashing == true ) {
					continue;
				}

				//create the minipose to use for early scoring if it does not exist already
				if ( minipose->size() == 0 ) {
					working_pose_->append_residue_by_jump(*ligresOP, 1);

					for ( core::Size resi_pos = 1; resi_pos < working_pose_->size(); ++resi_pos ) {
						//code breaks if unmatched disulfide  bonds form, just place all  residues that can have the disulfide type
						if ( working_pose_->residue(resi_pos).has_variant_type(core::chemical::DISULFIDE) ) {
							continue;
						}

						//code to try to make minipose even smaller to the point of only near residues
						//only consider adding residues that are within a distance equal to the sum of the nbr radius of the ligand (located at index pose.size) and residue being investigated
						if ( working_pose_->residue(working_pose_->size()).nbr_atom_xyz().distance(working_pose_->residue(resi_pos).nbr_atom_xyz()) < (working_pose_->residue(working_pose_->size()).nbr_radius() + working_pose_->residue(resi_pos).nbr_radius()) ) {
							//append residue to minipose
							minipose->append_residue_by_jump(working_pose_->residue(resi_pos), 1);
							ms_tr << resi_pos << ", ";
						}
					}

					//append ligand to minipose
					minipose->append_residue_by_jump(working_pose_->residue(working_pose_->size()), 1);
					ms_tr << "Made minipose of size " << minipose->size() << std::endl;

					//hard wipe minipose and then move to next placement if minipose only has the ligand in it
					if ( minipose->size() == 1 ) {
						core::pose::PoseOP filler(new pose::Pose);
						minipose = filler;
						//wipe ligand from working pose since we are not investigating this placement
						working_pose_->delete_residue_slow(working_pose_->size());
						continue;
					}

					//dump the minipose for debugging
					core::io::pdb::dump_pdb( *minipose, "minipose.pdb");
					//delete ligand so we can reuse minipose and wipe from working pose
					minipose->delete_residue_slow(minipose->size());
					working_pose_->delete_residue_slow(working_pose_->size());
				}

				//append ligand to minipose for early scoring
				minipose->append_residue_by_jump(*ligresOP, 1);

				fa_rep_fxn_->score(*minipose);

				//high fa_rep means clashing, want low fa_rep
				core::Real fa_rep = minipose->energies().residue_total_energies(minipose->size())[core::scoring::fa_rep];

				//check if fa_rep is good
				//positive score is bad
				//best scores are negative and closest to 0
				if ( fa_rep > fa_rep_cutoff ) {
					minipose->delete_residue_slow(minipose->size());
					++clashing_counter;
					continue;
				}

				fa_atr_fxn_->score(*minipose);

				core::Real fa_atr = minipose->energies().residue_total_energies(minipose->size())[core::scoring::fa_atr];

				//run fa_atr check
				//don't keep if fa_atr is greater than cutoff (notable values are -5 and -3.77)
				if ( fa_atr > fa_atr_cutoff ) {
					minipose->delete_residue_slow(minipose->size());
					++clashing_counter;
					continue;
				}
				++passed_placement_counter;

				//append ligand to working pose
				//using version of append to put it in a new chain
				working_pose_->append_residue_by_jump(*ligresOP, working_pose_->size(), "", "", true);

				//create constraints for ligand wiggling to add to the pose
				constraints::ConstraintSetOP sc_cst_set( new constraints::ConstraintSet() );

				core::scoring::func::FuncOP fx1( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
				sc_cst_set->add_constraint( core::scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_1, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position_ ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_1 ), fx1 ) ) );

				core::scoring::func::FuncOP fx2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
				sc_cst_set->add_constraint( core::scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_2, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position_ ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_2 ), fx2 ) ) );

				core::scoring::func::FuncOP fx3( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
				sc_cst_set->add_constraint( core::scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_3, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position_ ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_3 ), fx3 ) ) );

				working_pose_->constraint_set(sc_cst_set);

				//get free energy of pose with placed ligand
				std::map< std::string, core::Real > interface_mapX = protocols::ligand_docking::get_interface_deltas('2', *working_pose_, whole_score_fxn_no_constraints_, "");

				core::Real delta_score = interface_mapX["interface_delta_2"];

				whole_score_fxn_->score(*working_pose_);
				fa_rep = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_rep];
				fa_atr = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_atr];
				core::Real sc_constraint_check( working_pose_->energies().total_energies()[ coordinate_constraint ] );

				//attempt to use movers to optimize placement a little more
				ms_tr << "Pre-move delta score = " << delta_score << ", fa_atr = " << fa_atr << ", fa_rep = " << fa_rep << ", coordinate_constraint = " << sc_constraint_check << std::endl;

				//make ligandareaop for use with the highresdocker
				protocols::ligand_docking::LigandAreaOP sc_ligand_area(new protocols::ligand_docking::LigandArea());
				//adjust values for the LigandArea object (doesn't look like there is a constructor for it, but it has free access to variables)
				//using values from integration test for 7cpa ligand docking from xml file when applicable

				sc_ligand_area->chain_ = '^';
				sc_ligand_area->cutoff_ = 1;
				sc_ligand_area->add_nbr_radius_ = true;
				sc_ligand_area->all_atom_mode_ = true;
				sc_ligand_area->minimize_ligand_ = 1;

				core::pose::PDBInfoCOP pdb_info( working_pose_->pdb_info() );

				//add ligand to pose residue type set if it is not already in the set
				//program doesn't work if the ligand isn't added to the pose residue type sets
				if ( ligand_added == false ) {
					ligand_added = true;

					//code to add type set of imported ligand into pose
					core::chemical::PoseResidueTypeSetOP rts( working_pose_->conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );

					//check if ligand is already in rts. if it is not, then add it. otherwise continue (if it is already in, trying to add again breaks things)
					//core::chemical::MutableResidueTypeOP type_to_add( new core::chemical::MutableResidueType( pose->residue(pose->size()).type() ) );

					//make a ResidueTypeSetCOP that will be pulled from the PoseResidueTypeSetOP
					core::chemical::ResidueTypeSetCOP def_rts(rts->default_rts());

					//use the name_mapOP function to get a residue type pointer based on the name of the ligand (could be either a nullptr or a pointer to the type)
					//should be a nullptr if it isn't in the set
					core::chemical::ResidueTypeCOP lig_rt(def_rts->name_mapOP(working_pose_->residue(working_pose_->size()).name()));

					if ( lig_rt == nullptr ) {
						rts->add_base_residue_type(lig_mrt);
					}
					working_pose_->conformation().reset_residue_type_set_for_conf(rts);
				}

				//add ligand_area to an op of ligandarea
				utility::vector1<protocols::ligand_docking::LigandAreaOP> ligand_areas;
				ligand_areas.push_back(sc_ligand_area);

				//set up interfaces and movemaps for the highresdocker

				//make interfacebuilder
				protocols::ligand_docking::InterfaceBuilderOP sc_interface(new protocols::ligand_docking::InterfaceBuilder(ligand_areas));

				//make a blank interfacebuilder since the movemapbuilder needs 2
				protocols::ligand_docking::InterfaceBuilderOP blank_interface(new protocols::ligand_docking::InterfaceBuilder());


				//create movemapbuilderOP
				//we will minimize waters (boolean)
				protocols::ligand_docking::MoveMapBuilderOP my_movemapbuilder ( new protocols::ligand_docking::MoveMapBuilder(sc_interface, blank_interface, true));

				//use movemapbuilder to make highresdocker
				protocols::ligand_docking::HighResDockerOP my_HighResDocker( new protocols::ligand_docking::HighResDocker(3, 0, fa_atr_rep_fxn_, my_movemapbuilder) );

				//set highresdocker to not repack
				my_HighResDocker->set_allow_repacking(false);

				//apply with highresdocker
				my_HighResDocker->apply(*working_pose_);

				//score pose after using highresdocker
				std::map< std::string, core::Real > interface_mapX_postdock = protocols::ligand_docking::get_interface_deltas('2', *working_pose_, whole_score_fxn_no_constraints_, "");

				delta_score = interface_mapX_postdock["interface_delta_2"];
				//delta_score = delta_score * -1;

				whole_score_fxn_->score(*working_pose_);
				fa_rep = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_rep];
				fa_atr = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_atr];
				sc_constraint_check = working_pose_->energies().total_energies()[ coordinate_constraint ] ;

				ms_tr << "Post-dock delta score = " << delta_score << ", fa_atr = " << fa_atr << ", fa_rep = " << fa_rep << ", coordinate_constraint = " << sc_constraint_check << std::endl;

				//check if whole_score is within cutoff, kill if not
				//need to remove ligand from poses so that they can be recycled

				if ( delta_score > score_cutoff || fa_atr > fa_atr_cutoff || fa_rep > fa_rep_cutoff ) {
					minipose->delete_residue_slow(minipose->size());
					working_pose_->delete_residue_slow(working_pose_->size());
					++clashing_counter;
					continue;
				}


				//name the pdb  that could come from the pose
				//current naming convention
				std::string pdb_name = output_prefix + "_ResPos_" + std::to_string(working_position_) + "_ResID_" + discovery_position_residue + "_Trio" + std::to_string(i) + "_" + ligresOP->name() + "_motif_" + motifcop->remark() + "_rep_" + std::to_string(fa_rep) + "_atr_" + std::to_string(fa_atr) + "_delta_" + std::to_string(delta_score) + "_constr_" + std::to_string(sc_constraint_check) + ".pdb";


				//if we are keeping all placements that are better than a given ddg cutoff, it is better to not inflate the best_100_placements vector since we only use it for sorting
				//instead, just make the pdb and keep going
				if ( best_pdbs_to_keep == 0 ) {
					ms_tr << "Making pdb file for: " << pdb_name << std::endl;
					core::io::pdb::dump_pdb(*working_pose_, pdb_name);

					//directly remove residue from end of pose
					minipose->delete_residue_slow(minipose->size());
					working_pose_->delete_residue_slow(working_pose_->size());
					++passing_counter;

					//continue since we do not need to touch best_100_placements
					continue;
				}

				//if we are keeping a set number of placements, then we need to append to a list that will be sorted
				std::tuple<core::Real, core::pose::Pose, std::string> pose_tuple(delta_score, *working_pose_, pdb_name);

				//directly remove residue from end of pose
				minipose->delete_residue_slow(minipose->size());

				working_pose_->delete_residue_slow(working_pose_->size());

				//add pose tuple, heapify if we want to only keep the top X placements (determined by using the best_pdbs_to_keep != 0)
				//push pose into vector of passing placements
				if ( best_placements.size() < best_pdbs_to_keep ) {
					best_placements.push_back(pose_tuple);
					++passing_counter;
					if ( best_placements.size() == best_pdbs_to_keep ) {
						std::make_heap(best_placements.begin(), best_placements.end(), comparator_v);

					}
				} else if ( std::get<0>(best_placements.front()) > std::get<0>(pose_tuple) ) {
					std::pop_heap (best_placements.begin(),best_placements.end(), comparator_v); best_placements.pop_back();
					best_placements.push_back(pose_tuple); std::push_heap (best_placements.begin(),best_placements.end(), comparator_v);
				}

			}

			ms_tr << "# passing cases for this trio = " << passing_counter << std::endl;
			ms_tr << "# clashing cases for this trio = " << clashing_counter << std::endl;
			ms_tr << "# post motif filter cases for this trio = " << pose_atom_check_counter << std::endl;

			ligand_passing_counter += passing_counter;
			ligand_clashing_counter += clashing_counter;

			ms_tr << "Total cases for trio: " << passing_counter + clashing_counter << std::endl;

		}//end of trio loop

		ms_tr << "Done iterating all trios, moving to next ligand" << std::endl;
		ms_tr << "Total passing attempts for ligand is " << ligand_passing_counter << std::endl;
		ms_tr << "Total clashing attempts for ligand is " << ligand_clashing_counter << std::endl;


		//assign the new cutoff once we hit 100 entries
		if ( best_placements.size() >= best_pdbs_to_keep && best_pdbs_to_keep != 0 ) {
			//score_cutoff = std::get<0>(best_100_placements[best_100_placements.size() - 1]);
			score_cutoff = std::get<0>(best_placements[0]);
			ms_tr << "New score cutoff is: " << score_cutoff << std::endl;
		}
	}
	//sort passing_placements if keeping only a specified amount
	if ( best_pdbs_to_keep != 0 ) {
		best_placements = mergeSort(best_placements);
	}

	//create pdbs of best scoring poses for pdb

	for ( core::Size best_pose_runner = 0; best_pose_runner < best_placements.size(); best_pose_runner++ ) {
		ms_tr << "Making pdb file for: " << std::get<2>(best_placements[best_pose_runner]) << std::endl;
		core::io::pdb::dump_pdb( ( std::get<1>(best_placements[best_pose_runner])), std::get<2>(best_placements[best_pose_runner]));
	}

	ms_tr << "Number of placements that passed all filters: " << passed_placement_counter << std::endl;

	return 1;
}

//function to get a sub-library of motifs from the main library, based on the residue being used (only get for select residue)
//function may have use outside discover, so public can use it
protocols::motifs::MotifCOPs LigandDiscoverySearch::get_motif_sublibrary_by_aa(std::string residue_name)
{
	static basic::Tracer ms_tr( "LigandDiscoverySearch_get_motif_sublibrary_by_aa", basic::t_info );

	//create temporary motifcops object to hold
	protocols::motifs::MotifCOPs motif_holder;

	for ( auto motifcop : motif_library_ ) {
		//collect residue name from motifcop
		std::string motif_residue_name(motifcop->restype_name1());

		//if residue name matches that of the residue at the working position, add it to the motif_holder
		if ( motif_residue_name == residue_name ) {
			motif_holder.push_back(motifcop);
		}
	}

	//assign motif_holder contents to motif_library_for_select_residue_
	motif_library_for_select_residue_ = motif_holder;

	ms_tr << "Created motif sub-library for residue " << residue_name << " with " << motif_library_for_select_residue_.size() << " motifs in it." << std::endl;

	return motif_holder;
}

//create protein_representation_matrix_
//uses working_pose to make the matrix
void LigandDiscoverySearch::create_protein_representation_matrix(core::Size & x_shift, core::Size & y_shift, core::Size & z_shift, int & x_bound_int, int & y_bound_int, int & z_bound_int)
{

	//create tracer to identify points of the run
	static basic::Tracer ms_tr( "LigandDiscoverySearch_create_protein_matrix", basic::t_info );

	//run through all atoms to derive a range of dimensions to contain the protein in a 3D  space
	//since we can't have negative indices, we need to normalize the coordinate values so that everything is positive
	//derive constant values based  on the most negative values in each dimension, and then add that constant to all coordinates

	int smallest_x = 1;
	int smallest_y = 1;
	int smallest_z = 1;

	int largest_x = 1;
	int largest_y = 1;
	int largest_z = 1;

	//create a list of coordinates of each atom to hold and work with to fill the protein_representation_matrix
	//can't seem to make a vector of xyzVector objects, so will need to just make a custome 2D vector  to  hold the data
	utility::vector1<numeric::xyzVector<int>> atom_coordinates;

	//determine largest and smallest x,y,z  values to determine dimensions of matrix
	for ( core::Size res_num = 1; res_num <= working_pose_->size(); ++res_num ) {
		for ( core::Size atom_num = 1; atom_num <= working_pose_->residue(res_num).natoms(); ++atom_num ) {
			//get the x,y,z data of the atom
			numeric::xyzVector<int> atom_xyz = working_pose_->residue(res_num).xyz(atom_num);
			//convert the coordinates to integers
			atom_xyz.x() = static_cast<int>(atom_xyz.x());
			atom_xyz.y() = static_cast<int>(atom_xyz.y());
			atom_xyz.z() = static_cast<int>(atom_xyz.z());

			//determine if any of the values  are the smallest
			if ( smallest_x > atom_xyz.x() ) {
				smallest_x = atom_xyz.x();
			}
			if ( smallest_y > atom_xyz.y() ) {
				smallest_y = atom_xyz.y();
			}
			if ( smallest_z > atom_xyz.z() ) {
				smallest_z = atom_xyz.z();
			}

			//determine if any  are the largest
			if ( largest_x < atom_xyz.x() ) {
				largest_x = atom_xyz.x();
			}
			if ( largest_y < atom_xyz.y() ) {
				largest_y = atom_xyz.y();
			}
			if ( largest_z < atom_xyz.z() ) {
				largest_z = atom_xyz.z();
			}

			atom_coordinates.push_back(atom_xyz);

		}
	}

	//take negative values of the smallest values and then add 1 to derive the constants
	x_shift = (smallest_x * -1) + 1;
	y_shift = (smallest_y * -1) + 1;
	z_shift = (smallest_z * -1) + 1;

	//apply shift values to largest to get boundaries
	core::Size x_bound  = x_shift + largest_x;
	core::Size y_bound  = y_shift + largest_y;
	core::Size z_bound  = z_shift + largest_z;

	x_bound_int = x_bound;
	y_bound_int = y_bound;
	z_bound_int = z_bound;

	//apply constants to all coordinates
	for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates.size(); ++xyzVec ) {
		atom_coordinates[xyzVec].x() += x_shift;
		atom_coordinates[xyzVec].y() += y_shift;
		atom_coordinates[xyzVec].z() += z_shift;
	}

	//create 3D matrix to roughly represent 3D coordinate space of protein
	ms_tr << "Creating protein clash coordinate matrix. Dimensions of matrix are " << x_bound << "," << y_bound << "," << z_bound << std::endl;

	for ( core::Size x = 1; x <= x_bound; ++x ) {
		//make a 2D  matrix
		utility::vector1<utility::vector1<bool>> sub_matrix;

		for ( core::Size y = 1; y <= y_bound; ++y ) {

			//make a 1D matrix, seed with false values
			utility::vector1<bool> sub_sub_matrix(z_bound,  false);
			//push 1D  matrix into 2D
			sub_matrix.push_back(sub_sub_matrix);

		}
		//push a 2D  matrix into the 3D matrix
		protein_representation_matrix_.push_back(sub_matrix);
	}

	//seed the matrix with approximate coordinates of each atom
	//approximated by always rounding down via int casting
	for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates.size(); ++xyzVec ) {
		protein_representation_matrix_[static_cast<int>(atom_coordinates[xyzVec].x())][static_cast<int>(atom_coordinates[xyzVec].y())][static_cast<int>(atom_coordinates[xyzVec].z())] = true;
	}
}

//function to run a clash check of the placed ligand against the protein representation matrix
//ligand needs to be placed with a motif (at the very least needs coordinates)
bool LigandDiscoverySearch::ligand_clash_check(core::conformation::ResidueOP ligresOP, core::Size x_shift, core::Size y_shift, core::Size z_shift, int x_bound_int, int y_bound_int, int z_bound_int)
{
	//bool to track if ligand clashes
	bool has_clashing = false;

	//iterate through all atoms in the placed ligand, and determine if there is clashing
	//hold the number of times that there is clashing. If ligand clashes (ligand and protein atom in same cell), we can kill the attempt

	core::Size num_atoms_in_ligand = ligresOP->natoms();

	for ( core::Size residue_atom_iterator = 1; residue_atom_iterator <= num_atoms_in_ligand; ++residue_atom_iterator ) {
		//convert coordinates of current atom into format that can be read into protein placement matrix, determine if there is clashing
		numeric::xyzVector<int> test_atom_xyz = ligresOP->xyz(residue_atom_iterator);

		test_atom_xyz.x() = static_cast<int>(test_atom_xyz.x());
		test_atom_xyz.y() = static_cast<int>(test_atom_xyz.y());
		test_atom_xyz.z() = static_cast<int>(test_atom_xyz.z());

		//apply x,y,z shift to the coordinates
		test_atom_xyz.x() += x_shift;
		test_atom_xyz.y() += y_shift;
		test_atom_xyz.z() += z_shift;

		//handle case of ligand atom existing beyond protein matrix (won't be a clash anyway)
		if ( test_atom_xyz.x() < 1 || test_atom_xyz.x() > x_bound_int ) {
			continue;
		}
		if ( test_atom_xyz.y() < 1 || test_atom_xyz.y() > y_bound_int ) {
			continue;
		}
		if ( test_atom_xyz.z() < 1 || test_atom_xyz.z() > z_bound_int ) {
			continue;
		}

		//probe corresponding index of protein_representation_matrix
		if ( protein_representation_matrix_[test_atom_xyz.x()][test_atom_xyz.y()][test_atom_xyz.z()] ) {
			//clash case
			//since we hit a clash, no need to look at other ligand atoms
			has_clashing = true;
			break;
		}
	}

	return has_clashing;
}
