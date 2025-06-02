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
#include <core/scoring/etable/Etable.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>

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
#include <map>
#include <queue>
#include <functional>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
#include <basic/options/keys/protein_grid.OptionKeys.gen.hh>

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

#include <protocols/motifs/IdentifyLigandMotifs.hh>

#include <utility/Binary_Util.hh>

#include <protocols/protein_grid/ProteinGrid.hh>

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

//declare global tracer for class
static basic::Tracer ms_tr( "protocols.motifs.LigandDiscoverySearch" );

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
LigandDiscoverySearch::LigandDiscoverySearch(core::pose::PoseOP pose_from_PDB, protocols::motifs::MotifCOPs motif_library, utility::vector1<core::conformation::ResidueOP> all_residues, utility::vector1<core::Size> working_position)
{
	working_pose_ = pose_from_PDB;
	motif_library_ = motif_library;
	all_residues_ = all_residues;
	working_positions_ = working_position;

	//seed cutoff values
	seed_cutoff_values();

	//set up score functions
	setup_score_functions();
}

// @brief this function is to be called by the constructor(s) to seed initial values to cutoffs that are used for scoring/evaluating metrics of placed ligands in discover() and the functions it calls
void LigandDiscoverySearch::seed_cutoff_values()
{
	//score function cutoffs
	ms_tr.Debug << "Using fa_rep cutoff of: "  << option[ OptionKeys::motifs::fa_rep_cutoff ] << std::endl;
	fa_rep_cutoff_ = option[ OptionKeys::motifs::fa_rep_cutoff ];

	ms_tr.Debug << "Using fa_atr cutoff of: "  << option[ OptionKeys::motifs::fa_atr_cutoff ] << std::endl;
	fa_atr_cutoff_ = option[ OptionKeys::motifs::fa_atr_cutoff ];

	ms_tr.Debug << "Using fa_atr_rep cutoff of: "  << option[ OptionKeys::motifs::fa_atr_rep_cutoff ] << std::endl;
	fa_atr_rep_cutoff_ = option[ OptionKeys::motifs::fa_atr_rep_cutoff ];

	ms_tr.Debug << "Using whole score function cutoff of: "  << option[ OptionKeys::motifs::ligand_wts_fxn_cutoff ] << std::endl;
	whole_fxn_cutoff_ = option[ OptionKeys::motifs::ligand_wts_fxn_cutoff ];

	ms_tr.Debug << "Using ddg cutoff of: "  << option[ OptionKeys::motifs::ddg_cutoff ] << std::endl;
	ddg_cutoff_ = option[ OptionKeys::motifs::ddg_cutoff ];

	//placement motifs cutoffs
	ms_tr.Debug << "Using minimum motifs-like interactions cutoff of: "  << option[ OptionKeys::motifs::minimum_motifs_formed_cutoff] << std::endl;
	min_motifs_cutoff_ = option[ OptionKeys::motifs::minimum_motifs_formed_cutoff];

	ms_tr.Debug << "Using minimum motifs-like interactions on significant residues cutoff of: "  << option[ OptionKeys::motifs::minimum_significant_motifs_formed_cutoff] << std::endl;
	min_sig_motifs_cutoff_ = option[ OptionKeys::motifs::minimum_significant_motifs_formed_cutoff];

	ms_tr.Debug << "Using minimum real motifs interactions ratio cutoff of: "  << option[ OptionKeys::motifs::minimum_ratio_of_real_motifs_from_ligand] << std::endl;
	real_motif_ratio_cutoff_ = option[ OptionKeys::motifs::minimum_ratio_of_real_motifs_from_ligand];
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

//function to define the vector of residue indices that we will use for applying motifs for ligand placements
void LigandDiscoverySearch::set_working_positions(utility::vector1<core::Size> working_position)
{
	working_positions_ = working_position;
}

//overload of set_working_position to take in a single core::Size that is not in a vector, in case the user is only interesting in one position
//This seems like a handy way to be more flexible and allow working on a single position that doesn't have to be put in a vector first
void LigandDiscoverySearch::set_working_positions(core::Size working_position)
{
	//wipe previous contents of working_positions_
	working_positions_.clear();

	//push back the passed working_position
	working_positions_.push_back(working_position);
}

//append additional positions to working_positions_ to the end of the vector
void LigandDiscoverySearch::add_working_positions(utility::vector1<core::Size> working_positions)
{

	for ( const auto & working_position : working_positions ) {
		//push back the passed working_position
		working_positions_.push_back(working_position);
	}
}

//append additional position to working_positions_ to the end of the vector
void LigandDiscoverySearch::add_working_positions(core::Size working_position)
{

	//push back the passed working_position
	working_positions_.push_back(working_position);

}

//return contents of working_positions_
//probably always safest to return as a vector, even if there could only be one entry
utility::vector1<core::Size> LigandDiscoverySearch::get_working_positions()
{
	return working_positions_;
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
	pose::PoseOP pose_pointer_helper(pose_helper.clone());
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

// @brief prepare score functions for usage in discovery function. Called within discover() and not in a constructor. This probably shouldn't be messed with, so it is kept private
void LigandDiscoverySearch::setup_score_functions()
{
	//initially seed score class score functions
	whole_score_fxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "ligand.wts" );
	fa_atr_fxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "ligand.wts" );
	fa_rep_fxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "ligand.wts" );
	fa_atr_rep_fxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "ligand.wts" );

	//for each weight in the whole score function, set the scoretype weight to 0
	//the purpose of this is to set up the non-weight parameters exactly the same way (and using the terms from the ligand.wts function)
	//once terms are blanked to 0, terms of interest will be re-weighted as seen below
	for ( auto my_scoretype : whole_score_fxn_->get_nonzero_weighted_scoretypes() ) {
		fa_rep_fxn_->set_weight(my_scoretype,0);
		fa_atr_fxn_->set_weight(my_scoretype,0);
		fa_atr_rep_fxn_->set_weight(my_scoretype,0);
	}

	//set weight values for atr and rep based on values in ligand.wts
	fa_rep_fxn_->set_weight(core::scoring::fa_rep, 0.8);

	fa_atr_fxn_->set_weight(core::scoring::fa_atr, 0.4);

	fa_atr_rep_fxn_->set_weight(core::scoring::fa_atr, 0.8);
	fa_atr_rep_fxn_->set_weight(core::scoring::fa_rep, 0.4);

	//seed working function to be either the whole or atr_rep function, which will be used on most heavy duty scoring cases
	//there are individual cases where the user can additionally select to get a score with a specific function that will be left as is, but this working_fxn will be used in all other cases
	if ( option[ OptionKeys::motifs::highresdock_with_whole_score_fxn ] ) {
		working_fxn_ = whole_score_fxn_;
	} else {
		working_fxn_ = fa_atr_rep_fxn_;
	}

	//set coordinate_constraint to zero for the whole function (creates strange ddg scores otherwise)
	whole_score_fxn_->set_weight(core::scoring::coordinate_constraint, 0);
}

// @brief function to push all adjacent atom indices if inputted ligand residue (pointer) into a 3D core::Size vector (atom_trios typedef)
LigandDiscoverySearch::atom_trios LigandDiscoverySearch::derive_adjacent_atoms_of_ligand(const core::conformation::ResidueOP ligresOP, const core::chemical::AtomTypeSetCOP atset)
{
	atom_trios ligand_atom_trios;

	//find all atom trios (that do not contain hydrogen) in the ligand
	for ( core::Size atom_i = 1; atom_i <= ligresOP->natoms(); ++atom_i ) {
		if ( ligresOP->atom_is_hydrogen(atom_i) ) { continue; }
		// This is a for loop to iterate over each atom's connected atoms:
		core::conformation::Residue::AtomIndices atom_i_connects(  ligresOP->bonded_neighbor( atom_i ) );
		for ( const auto & atom_j :  atom_i_connects ) {
			if ( ligresOP->atom_is_hydrogen(atom_j) ) { continue; }
			// This is the next for loop to find connects for the second atom, giving us the final atom number (atom k)
			core::conformation::Residue::AtomIndices atom_j_connects(  ligresOP->bonded_neighbor( atom_j ) );
			for ( const auto & atom_k :  atom_j_connects ) {
				if ( ligresOP->atom_is_hydrogen(atom_k) ) { continue; }
				chemical::AtomType atom_i_type(ligresOP->atom_type(atom_i));
				if ( atom_i != atom_k ) {

					//make the 3 atom vector
					utility::vector1< utility::vector1< core::Size > > cur_motif_indices;

					utility::vector1< core::Size > atom_i_vector;
					atom_i_vector.push_back( atom_i );
					atom_i_vector.push_back( atset->atom_type_index( ligresOP->atom_type(atom_i).atom_type_name() ) );

					utility::vector1< core::Size > atom_j_vector;
					atom_j_vector.push_back( atom_j );
					atom_j_vector.push_back( atset->atom_type_index( ligresOP->atom_type(atom_j).atom_type_name() ) );

					utility::vector1< core::Size > atom_k_vector;
					atom_k_vector.push_back( atom_k );
					atom_k_vector.push_back( atset->atom_type_index( ligresOP->atom_type(atom_k).atom_type_name() ) );

					cur_motif_indices.push_back( atom_i_vector);
					cur_motif_indices.push_back( atom_j_vector);
					cur_motif_indices.push_back( atom_k_vector);

					ligand_atom_trios.push_back(cur_motif_indices);
				}
			}
		}
	}

	return ligand_atom_trios;

}


// @brief function used to make a minipose (focused pose around placed ligand to get quicker scoring of metrics like fa_atr and fa_rep)
//if returns true, a minipose was successfully made; if returns false, minipose is still empty because no other residues were recruited to it
bool LigandDiscoverySearch::make_minipose(core::pose::PoseOP & minipose, const core::conformation::ResidueOP ligresOP)
{
	working_pose_->append_residue_by_jump(*ligresOP, 1);

	ms_tr.Debug << "Builting Minipose made of residue indices: ";

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

			ms_tr.Debug << resi_pos << ", ";

		}
	}

	ms_tr.Debug << std::endl;

	//append ligand to minipose
	minipose->append_residue_by_jump(working_pose_->residue(working_pose_->size()), 1);

	ms_tr.Debug << "Made minipose of size " << minipose->size() << std::endl;


	//hard wipe minipose and then move to next placement if minipose only has the ligand in it
	if ( minipose->size() == 1 ) {
		core::pose::PoseOP filler(new pose::Pose);
		minipose = filler;
		//wipe ligand from working pose since we are not investigating this placement
		working_pose_->delete_residue_slow(working_pose_->size());
		//return false so discover() knows to continue to the next placement, becase the current is bad
		return false;
	}

	//dump the minipose for debugging
	core::io::pdb::dump_pdb( *minipose, "minipose.pdb");
	//delete ligand so we can reuse minipose and wipe from working pose
	minipose->delete_residue_slow(minipose->size());
	working_pose_->delete_residue_slow(working_pose_->size());
	return true;
}

// @brief function to score the minipose iteratively on fa_rep, fa_atr, and fa_atr+fa_rep
//this is performed iteratively in order to more quickly filter on select criteria with more and quicker filtering happening at earlier steps
//returns a boolean based on whether the minipose score was good enough at all steps or not (false means it failed, and the placement will be killed)
bool LigandDiscoverySearch::score_minipose(const core::pose::PoseOP & minipose, core::Real & fa_rep, core::Real & fa_atr, core::Real & fa_atr_rep_score_before)
{
	fa_rep_fxn_->score(*minipose);

	//high fa_rep means clashing, want low fa_rep
	fa_rep = minipose->energies().residue_total_energies(minipose->size())[core::scoring::fa_rep];

	ms_tr.Debug << "fa_rep = " << fa_rep << std::endl;

	//check if fa_rep is good
	//positive score is bad
	//best scores are negative and closest to 0
	if ( fa_rep > fa_rep_cutoff_ ) {
		return false;
	}

	fa_atr_fxn_->score(*minipose);

	fa_atr = minipose->energies().residue_total_energies(minipose->size())[core::scoring::fa_atr];

	ms_tr.Debug << "fa_atr = " << fa_rep << std::endl;

	//run fa_atr check
	//don't keep if fa_atr is greater than cutoff
	if ( fa_atr > fa_atr_cutoff_ ) {
		return false;
	}

	//score whole minipose with atr_rep function
	fa_atr_rep_score_before = fa_atr_rep_fxn_->score(*minipose);

	ms_tr.Debug << "fa_atr_rep = " << fa_rep << std::endl;

	//check if worse than cutoff
	if ( fa_atr_rep_score_before > fa_atr_rep_cutoff_ ) {
		return false;
	}

	return true;
}

// @brief create a constraint set on the 3 ligand motif atoms (the last residue in working_pose_) to reduce their movement before applying a highresdock
void LigandDiscoverySearch::add_constraints_to_working_pose(const core::Size trip_atom_1, const core::Size trip_atom_2, const core::Size trip_atom_3, const core::Size working_position, const core::conformation::ResidueOP ligresOP)
{
	//create constraints for ligand wiggling to add to the pose
	constraints::ConstraintSetOP sc_cst_set( new constraints::ConstraintSet() );

	core::scoring::func::FuncOP fx1( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	sc_cst_set->add_constraint( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_1, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_1 ), fx1 ) );

	core::scoring::func::FuncOP fx2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	sc_cst_set->add_constraint(  utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_2, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_2 ), fx2 ) );

	core::scoring::func::FuncOP fx3( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	sc_cst_set->add_constraint( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_3, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_3 ), fx3 ) );

	working_pose_->constraint_set(sc_cst_set);
}

// @brief this function takes in a selected scorefunctionOP and gets the ddg of the selected poseOP (ideally with a placed ligand), and returns the ddg
core::Real LigandDiscoverySearch::get_pose_ddg(core::scoring::ScoreFunctionOP score_fxn, core::pose::PoseOP & my_pose)
{
	//use the passed score function to get interface deltas on the working_pose
	std::map< std::string, core::Real > interface_mapX_postdock = protocols::ligand_docking::get_interface_deltas('2', *my_pose, score_fxn, "");

	//return the ddg, which is interface_delta_2 in the interface map
	return interface_mapX_postdock["interface_delta_2"];
}

// @brief This function adds a ligand mutableresiduetypeOP to the residue type set for working_pose and original_pose. This occurs on the first instance of this ligand passing enough filters to be used in highresdock
// the ligand is added to both so that the ligand can be a part of the original_pose for future iterations of placement attempts for the ligand, as well as the working_pose for immediate use
void LigandDiscoverySearch::add_ligand_to_pose_residuetypeset(const core::chemical::MutableResidueTypeOP lig_mrt)
{
	//code to add type set of imported ligand into pose
	core::chemical::PoseResidueTypeSetOP rts( working_pose_->conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t ) );

	//keeping this commented just in case, but now removing the pull of the default rts, since I now know you can just pull the ligand residue type straight from the poseRTS
	//make a ResidueTypeSetCOP that will be pulled from the PoseResidueTypeSetOP
	//core::chemical::ResidueTypeSetCOP def_rts(rts->default_rts());

	//use the name_mapOP function to get a residue type pointer based on the name of the ligand (could be either a nullptr or a pointer to the type)
	//should be a nullptr if it isn't in the set
	//core::chemical::ResidueTypeCOP lig_rt(def_rts->name_mapOP(working_pose_->residue(working_pose_->size()).name()));
	core::chemical::ResidueTypeCOP lig_rt(rts->name_mapOP(working_pose_->residue(working_pose_->size()).name()));

	//add the ligand to the residue type set if lig_rt is a nullpointer (which it should be)
	if ( lig_rt == nullptr ) {
		rts->add_base_residue_type(lig_mrt);
		//reset the residue type sets for the working_pose_ and the original_pose_ so that they are updated to have the new ligand
		working_pose_->conformation().reset_residue_type_set_for_conf(rts);
		original_pose_.conformation().reset_residue_type_set_for_conf(rts);
	}
}

// @brief This function creates a HighResDockOP object for using highresdock in discover() to try to optimize the ligand placement. The function takes in a score function OP to use in the HRD
protocols::ligand_docking::HighResDockerOP LigandDiscoverySearch::make_HighResDockOP_for_discovery(const core::scoring::ScoreFunctionOP my_fxn)
{
	//make ligandareaop for use with the highresdocker
	protocols::ligand_docking::LigandAreaOP sc_ligand_area(new protocols::ligand_docking::LigandArea());
	//adjust values for the LigandArea object (doesn't look like there is a constructor for it, but it has free access to variables)
	//using values from integration test for 7cpa ligand docking from xml file when applicable

	sc_ligand_area->chain_ = '^';
	sc_ligand_area->cutoff_ = 1;
	sc_ligand_area->add_nbr_radius_ = true;
	sc_ligand_area->all_atom_mode_ = true;
	sc_ligand_area->minimize_ligand_ = 1;

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

	//use movemapbuilder to make highresdocker; make 1 HRD for each score function
	//first 2 values correspond to the number of mover cycles and to repack every Nth cycle
	//we will move 3 times and do not want to repack (residues seem to move to unwanted positions)
	protocols::ligand_docking::HighResDockerOP my_HighResDocker( new protocols::ligand_docking::HighResDocker(3, 0, my_fxn, my_movemapbuilder) );

	//set highresdockers to not repack
	my_HighResDocker->set_allow_repacking(false);

	return my_HighResDocker;
}

// @brief This function uses a HighResDockOP object for using highresdock in discover() to try to optimize the ligand placement
void LigandDiscoverySearch::run_HighResDock_on_working_pose(const protocols::ligand_docking::HighResDockerOP my_HighResDocker)
{
	my_HighResDocker->apply(*working_pose_);
}

//main function to run ligand discovery operations
//needs to have values set for working_pose_, motif_library_, and all_residues_
//parameter is a string to be a prefix name to use for outputted file names
core::Size LigandDiscoverySearch::discover(std::string output_prefix)
{
	//create identifyligandmotifs object for use with deriving motifs from placed ligands
	IdentifyLigandMotifs ilm;

	// Make an atomtypeset to get atomtype integers for use in derive_adjacent_atoms_of_ligand
	core::chemical::AtomTypeSetCOP atset = core::chemical::ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );

	//declare the highresdockerop
	protocols::ligand_docking::HighResDockerOP my_HighResDocker;

	//create the highresdocker object using either the working function
	my_HighResDocker = make_HighResDockOP_for_discovery(working_fxn_);

	//iterate over all indices in working_positions_
	//if the size of working_positions is 0, return -1 because we want at least 1 index to work with
	if ( working_positions_.size() == 0 ) {

		ms_tr.Warning << "Length of working_positions_ vector is 0 and we have no positions to work with, killing the attempt now." << std::endl;

		return -1;
	}

	//run discovery over each position in working_positions_
	for ( const auto & working_position : working_positions_ ) {

		ms_tr << "Current anchor residue position: " << working_position << std::endl;


		//declare empty discovery_position_residue name that will be properly written in initialization
		std::string discovery_position_residue = "";

		//determine whether to kill due to bad initialization
		bool kill_bad_init = false;

		//kill if working position is out of bounds
		if ( working_position < 1 || working_position > working_pose_->size() ) {

			ms_tr.Warning << "Working position of " << working_position << " is invalid as it is not a valid index to access a residue in your pose of size " << working_pose_->size() << std::endl;

			kill_bad_init = true;
		} else {
			//inputs are good if you make it here

			//set discovery_position_residue
			discovery_position_residue = working_pose_->residue(working_position).name3();

			//get motif sublibrary
			motif_library_for_select_residue_ = get_motif_sublibrary_by_aa(discovery_position_residue);

			if ( motif_library_for_select_residue_.size() == 0 ) {
				ms_tr.Warning << "We have no motifs to work with here. Exiting function." << std::endl;
				kill_bad_init = true;
			}
		}

		if ( kill_bad_init == true ) {

			ms_tr.Warning << "We have at least 1 bad initial input, killing the attempt now." << std::endl;

			return -1;
		}

		//derive motif_library_for_select_residue_ from motif_library and residue in working_pose_ and index working_position

		motif_library_for_select_residue_ = get_motif_sublibrary_by_aa(working_pose_->residue(working_position).name3());

		//hash the motif library into a map of smaller motif lists based on the residue involved and 6 atom names/types
		//define new map
		std::map<protocols::motifs::motif_atoms,protocols::motifs::MotifCOPs> mymap;
		//only hash if we use the flag check_if_ligand_motifs_match_real, since this is currently the only thing we use it for
		if ( option[ OptionKeys::motifs::check_if_ligand_motifs_match_real] ) {
			protocols::motifs::hash_motif_library_into_map(motif_library_,mymap);
		}

		//declare the clash pose grid

		//first, declare the resolution factor for the matrix
		core::Real resolution_increase_factor = option[ OptionKeys::motifs::resolution_scale_factor_float ];

		//next, identify if we want to define a sub-area (binding pocket) to explore further
		//we will set up the ProteinGrid depending on whether we have the sub area or not
		//use if we used binding_pocket_dimensions_sf and binding_pocket_center_sf flags to define the sub area dimensions and center
		if ( option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].user() && option[ OptionKeys::motifs::binding_pocket_center_sf ].user() ) {
			//make sure dimensions vector is valid
			//can not use if fewer than 3 values (do not use), warn if there are more than 3
			if ( option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size()<3 ) {
				ms_tr.Warning << "Input vector for binding_pocket_dimensions_sf is too small with only " << option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size() << " entries. We will not use a sub-area/binding pocket in analyses." << std::endl;
			} else {
				if ( option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size()>3 ) {
					ms_tr.Warning << "Input vector for binding_pocket_dimensions_sf has more than 3 entries with having " << option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size() << " entries. We will only use the first 3 values." << std::endl;
				}

				//create vector to hold the dimensions
				//store as a size
				utility::vector1<core::Size> binding_pocket_dimensions;

				//set the first 3 values of the binding_pocket_dimensions_sf vector to binding_pocket_dimensions (and also apply the resolution factor)
				binding_pocket_dimensions = option[ OptionKeys::motifs::binding_pocket_dimensions_sf ]();

				//now derive the binding pocket center, and ensure the input is valid
				numeric::xyzVector<int> bp_xyz(0, 0, 0);
				utility::vector1<int> binding_pocket_center_xyz = option[ OptionKeys::motifs::binding_pocket_center_sf ]();

				//ensure that there are enough entries in the center vector
				if ( binding_pocket_center_xyz.size() <= 2 ) {
					//handle if too small
					ms_tr.Warning << "You have only inputted " << binding_pocket_center_xyz.size() << " coordinates for the binding pocket center xyz coordinates. We need exactly 3 to work. We can't do anything with this, and will skip this method and will not use a sub-area/binding pocket in analyses." << std::endl;
					//make the proteingrid that just wraps around the working_pose_
					clash_pose_grid_ = utility::pointer::make_shared<protocols::protein_grid::ProteinGrid>(working_pose_, resolution_increase_factor);
				} else {
					//throw potential warnings if there are more than 3 entries in the input, and then make the grid that also has a sub area
					for ( core::Size coordinate = 1; coordinate <= binding_pocket_center_xyz.size(); ++coordinate ) {
						//apply relevant xyz shift as long as the coordinate is between 1-3, otherwise throw a warning
						if ( coordinate > 3 ) {

							ms_tr.Warning << "You have inputted an excess value using the motifs::binding_pocket_center_sf flag. This is value #" << coordinate << " in the input vector with a value of: " << binding_pocket_center_xyz[coordinate] << ". Ignoring this bad input, and you should review the proper usage for this flag." << std::endl;

						} else {
							//translation of contents of the input vector1 into a xyzvector to use in the ProteinGrid
							if ( coordinate == 1 ) {
								bp_xyz.x() = binding_pocket_center_xyz[coordinate];
							}
							if ( coordinate == 2 ) {
								bp_xyz.y() = binding_pocket_center_xyz[coordinate];
							}
							if ( coordinate == 3 ) {
								bp_xyz.z() = binding_pocket_center_xyz[coordinate];
							}
						}
					}

					//make the proteingrid that wraps around the working_pose_ and declares a sub area
					clash_pose_grid_ = utility::pointer::make_shared<protocols::protein_grid::ProteinGrid>(working_pose_, resolution_increase_factor, bp_xyz, binding_pocket_dimensions);
				}
			}
		} else {
			//make the proteingrid that just wraps around the working_pose_
			clash_pose_grid_ = utility::pointer::make_shared<protocols::protein_grid::ProteinGrid>(working_pose_, resolution_increase_factor);
		}

		//now create the space fill matrix, first as a copy of the clash matrix
		sf_pose_grid_ = utility::pointer::make_shared<protocols::protein_grid::ProteinGrid>(*clash_pose_grid_);

		//run space fill on the space fill grid
		sf_pose_grid_->project_lj_radii();

		//create a copy of the space fill matrix that will be used for checking the space fill with the ligand. The member sf_pose_grid will be used as a base to copy back over after the space fill analysis, since that should be faster than reverting the grid using built-in functions (since we have a copy to use as a shortcut)
		protocols::protein_grid::ProteinGridOP working_sf_pose_grid = utility::pointer::make_shared<protocols::protein_grid::ProteinGrid>(*sf_pose_grid_);

		//derive the occupied ratios of the empty system before placing any ligands, in case we are interested in deriving differentials in space fill
		//index 1 corresponds to the occupied ratio for the whole system
		//index 2 corresponds to only the sub-area
		utility::vector1<core::Real>  occupied_ratios (2,0);
		occupied_ratios[1] = working_sf_pose_grid->get_grid_occupied_cell_ratio();
		occupied_ratios[2] = working_sf_pose_grid->get_sub_grid_occupied_cell_ratio();

		//export the empty space fill matrix if selected
		if ( option[ OptionKeys::motifs::output_space_fill_matrix_pdbs ] ) {
			//std::string pdb_name = output_prefix + "_ResPos_" + std::to_string(working_position) + "_ResID_" + discovery_position_residue + "_Trio" + std::to_string(i) + "_" + ligresOP->name() + "_motif_" + motifcop->remark() + "_rep_" + std::to_string(fa_rep) + "_atr_" + std::to_string(fa_atr) + "_delta_" + std::to_string(delta_score) + "_constr_" + std::to_string(sc_constraint_check) + ".pdb";
			std::string matrix_pdb_prefix = output_prefix + "_empty_";
			//call function to export, use "empty" as the string prefix
			working_sf_pose_grid->export_protein_matrix_to_pdb(matrix_pdb_prefix);
		}


		//create and load in the space fill cutoff scores (if there are any), default to zero; zero can be an acceptable inputted score
		//store as a vector
		//index 1 is the cutoff for the whole system, index 2 is for the sub area, index 3 is for the sub area differential cutoff
		utility::vector1<core::Real> score_cutoffs_sf (3,0);

		if ( option[ OptionKeys::motifs::space_fill_cutoff_score ].user() ) {
			score_cutoffs_sf[1] = option[ OptionKeys::motifs::space_fill_cutoff_score ];
		}
		if ( option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() ) {
			score_cutoffs_sf[2] = option[ OptionKeys::motifs::space_fill_cutoff_score_sub ];
		}
		if ( option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {
			score_cutoffs_sf[3] = option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ];
		}

		//for whole ligand.wts function, determine also if we want to use the function for scoring at all
		bool use_ligand_wts = option[ OptionKeys::motifs::score_with_ligand_wts_function ];
		bool use_atr_rep = option[ OptionKeys::motifs::post_highresdock_fa_atr_rep_score ];

		//create vector to hold the top X placements
		comparator comparator_v = comparator();
		//this gets used unless the value for the number of best placements to collect is 0 (in which all are collected)
		std::vector < std::tuple<core::Real, core::pose::Pose, std::string>> best_placements;

		//determine how many output files to keep
		//if the value is 0, all placements that pass all filters will be kept
		core::Size best_pdbs_to_keep = 0;
		if ( option[ OptionKeys::motifs::best_pdbs_to_keep ].user() ) {
			best_pdbs_to_keep = option[ OptionKeys::motifs::best_pdbs_to_keep ];
		}

		ms_tr << "Starting to iterate through all ligands" << std::endl;

		//create a clone of the working_pose_ that working_pose_ can be reset to after each placement attempt
		original_pose_ = *((*working_pose_).clone());

		//hold the number of placements that pass all filters and could enter the top 100 placements
		int passed_placement_counter = 0;

		//now we have the filtered motif library to work with, run through each  atom trio in each ligand and try to match it against all motifs for the residue
		for ( const auto & tracker : all_residues_ ) {

			//counter to count the placements for this ligand that get passed to a pdb (used in clean naming)
			core::Size unique_placement_counter = 0;

			bool ligand_added = false;

			//core::chemical::ResidueTypeCOP ligres(ref);
			//convert ligres to be a ResidueOP type
			core::conformation::ResidueOP ligresOP = tracker;

			ms_tr << "On ligand " << ligresOP->name() << std::endl;

			const core::Real lig_nbr_radius = ligresOP->nbr_radius();

			ms_tr.Debug << "NBR_RADIUS of ligand is: " << lig_nbr_radius << std::endl;

			int ligand_passing_counter  = 0;
			int ligand_clashing_counter = 0;

			//derive a mutable residue type, used in high res docker
			core::chemical::MutableResidueTypeOP lig_mrt( new core::chemical::MutableResidueType( ligresOP->type() ) );

			//3D vector to hold data of all connected atom trios
			//contents of vector are as follows:
			//ligand_atom_trios[trio identifier: 1-#trios][atom identifier within trio: 1-3][trio atom metadata]
			//metadata consists of index (position 1) and numerical atom_type_index (position 2)

			ms_tr.Trace << "Finding all atom trios for this ligand" << std::endl;
			atom_trios ligand_atom_trios = derive_adjacent_atoms_of_ligand(ligresOP, atset);

			//mini pose to represent a smaller region of the protein near the binding pocket for early energy calculations
			core::pose::PoseOP minipose(new pose::Pose);

			//now we have all trios, iterate through each trio and the motif sub-library

			ms_tr.Debug << "Number of unique atom trios in this ligand are: " << ligand_atom_trios.size() << std::endl;


			//iterate through each atom trio and try to use motifs to place the ligand
			for ( core::Size i = 1; i <= ligand_atom_trios.size(); ++i ) {


				ms_tr.Trace << "On trio # " << i << std::endl;


				ms_tr.Debug << "Trio is " << ligresOP->atom_name(ligand_atom_trios[i][1][1]) << " " << ligresOP->atom_name(ligand_atom_trios[i][2][1]) << " " << ligresOP->atom_name(ligand_atom_trios[i][3][1]) << std::endl;

				core::Size trip_atom_1(ligand_atom_trios[i][1][1]);
				core::Size trip_atom_2(ligand_atom_trios[i][2][1]);
				core::Size trip_atom_3(ligand_atom_trios[i][3][1]);

				int clashing_counter = 0;
				int passing_counter = 0;
				int pose_atom_check_counter = 0;

				//run through the motif sublibrary and attempt to pair the ligand trio to the target residue based  on existing motifs
				//iterate through motif sublibrary

				ms_tr.Debug << "Looking through all motifs against this trio" << std::endl;
				ms_tr.Debug << "#motifs = " << motif_library_for_select_residue_.size() << std::endl;

				int motif_counter = 0;

				for ( auto motifcop : motif_library_for_select_residue_ ) {
					++motif_counter;

					if ( motif_counter % 10000 == 0 ) {
						ms_tr.Trace << "On motif #" << motif_counter <<std::endl;
					}

					//compare the atoms on the ligand side of the motif to the atom trio; continue if not a match
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
					motifcop->place_residue( working_pose_->residue(working_position), *ligresOP, trip_atom_1, trip_atom_2, trip_atom_3 , true );

					//bool to determine if placed residue clashes against the backbone
					bool has_clashing = false;

					//check if the placement clashes
					has_clashing = clash_pose_grid_->placed_ligand_clash_analysis(ligresOP);

					//continue because we clash
					if ( has_clashing == true ) {
						continue;
					}

					//space fill analysis block
					//faster than using score functions
					//this is used to determine if enough of a defined binding pocked is filled by a placed ligand or not (ensures elimination of off-target placements)
					//only attempt this block if we request at least one of the score cutoffs for space filling

					if ( option[ OptionKeys::motifs::space_fill_cutoff_score ].user() || option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() || option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {

						//track the space fill score for the system and sub_area, hold as a 2d Real vector
						//index 1 is the score for the whole system, index 2 is the score for the sub area
						//default values are 0
						utility::vector1<core::Real> space_fill_scores (2,0);

						//run space fill analysis function
						working_sf_pose_grid->placed_ligand_space_fill_analysis(ligresOP);

						//extract the full and sub ratios
						space_fill_scores[1] = working_sf_pose_grid->get_grid_occupied_cell_ratio();
						space_fill_scores[2] = working_sf_pose_grid->get_sub_grid_occupied_cell_ratio();

						//derive a differential score for the sub area between the placed and empty system
						//only do anything if the user used the differencial cutoff score
						//differential defined as the placed fill ratio - empty fill ratio for the respective sub areas

						core::Real fill_differential = 0;

						if ( option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {
							fill_differential = space_fill_scores[2] - occupied_ratios[2];
						}


						//export the space fill matrix if selected
						if ( option[ OptionKeys::motifs::output_space_fill_matrix_pdbs ] ) {
							//std::string pdb_name = output_prefix + "_ResPos_" + std::to_string(working_position) + "_ResID_" + discovery_position_residue + "_Trio" + std::to_string(i) + "_" + ligresOP->name() + "_motif_" + motifcop->remark() + "_rep_" + std::to_string(fa_rep) + "_atr_" + std::to_string(fa_atr) + "_delta_" + std::to_string(delta_score) + "_constr_" + std::to_string(sc_constraint_check) + ".pdb";
							std::string matrix_pdb_prefix = output_prefix + "_ResPos_" + std::to_string(working_position) + "_ResID_" + discovery_position_residue + "_Trio" + std::to_string(i) + "_" + ligresOP->name() + "_motif_" + motifcop->remark();
							//call function to export, use "empty" as the string prefix
							working_sf_pose_grid->export_protein_matrix_to_pdb(matrix_pdb_prefix);
						}

						//at end before check, reset the working_sf_pose_grid to wipe the placed ligand
						working_sf_pose_grid = utility::pointer::make_shared<protocols::protein_grid::ProteinGrid>(*sf_pose_grid_);

						//run check for if the placement is passable based on score
						//either the whole system score or the sub area score need to pass (one can fail)
						if ( (space_fill_scores[1] >= score_cutoffs_sf[1] && option[ OptionKeys::motifs::space_fill_cutoff_score ].user()) || (space_fill_scores[2] >= score_cutoffs_sf[2] && option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user()) || (fill_differential >= score_cutoffs_sf[3] && option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user()) ) {
							//print out where any scores passed

							if ( space_fill_scores[1] >= score_cutoffs_sf[1] && option[ OptionKeys::motifs::space_fill_cutoff_score ].user() ) {
								ms_tr.Debug << space_fill_scores[1] << " : whole system, passed" << std::endl;
							} else if ( option[ OptionKeys::motifs::space_fill_cutoff_score ].user() ) {
								ms_tr.Debug << space_fill_scores[1] << " : whole system, failed" << std::endl;
							}

							if ( space_fill_scores[2] >= score_cutoffs_sf[2] && option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() ) {
								ms_tr.Debug << space_fill_scores[2] << " : sub-area, passed" << std::endl;
							} else if ( option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() ) {
								ms_tr.Debug << space_fill_scores[2] << " : sub-area, failed" << std::endl;
							}

							if ( fill_differential >= score_cutoffs_sf[3] && option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {
								ms_tr.Debug << fill_differential << " : sub-area differential, passed" << std::endl;
							} else if ( option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {
								ms_tr.Debug << fill_differential << " : sub-area differential, failed" << std::endl;
							}

						} else {
							//print for any failed scores

							if ( option[ OptionKeys::motifs::space_fill_cutoff_score ].user() ) {
								ms_tr.Debug << space_fill_scores[1] << " : whole system, failed and killing palcement attempt" << std::endl;
							}
							if ( option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() ) {
								ms_tr.Debug << space_fill_scores[2] << " : sub-area, failed and killing placement attempt" << std::endl;
							} else if ( option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {
								ms_tr.Debug << fill_differential << " : sub-area differential, failed" << std::endl;
							}

							continue;
						}
					}

					//create the minipose to use for early scoring (atr, rep, atrrep) if it does not exist already
					if ( minipose->size() == 0 ) {
						//try to make the minipose if it is currently empty
						//if the function returns false, that means that the minipose is bad and has no residues beyond the ligand; we want to continue if this is the case and the next placement will attempt
						if ( make_minipose(minipose, ligresOP) == false ) {
							continue;
						}
					}

					//append ligand to minipose for early scoring
					minipose->append_residue_by_jump(*ligresOP, 1);

					//set up values to be passed into score_minipose and used downstream
					core::Real fa_rep;
					core::Real fa_atr;
					core::Real fa_atr_rep_score_before;

					//score the minipose and set the value to a bool
					bool minipose_scoring = score_minipose(minipose,fa_rep,fa_atr,fa_atr_rep_score_before);



					//delete the last residue off minipose (the ligand), since we no longer need it and can recycle the minipose for further iterations
					minipose->delete_residue_slow(minipose->size());

					//if the value is true, the minipose passes initial scoring
					//if false, the placement can be killed
					if ( !minipose_scoring ) {
						++clashing_counter;
						continue;
					}

					++passed_placement_counter;

					//append ligand to working pose
					//using version of append to put it in a new chain
					working_pose_->append_residue_by_jump(*ligresOP, working_pose_->size(), "", "", true);

					//apply constraint set to working_pose_
					add_constraints_to_working_pose(trip_atom_1, trip_atom_2, trip_atom_3, working_position, ligresOP);

					//get free energy of pose with placed ligand before highresdock
					//declaration of variable
					core::Real delta_score = 10000;

					//use working function based to get ddg before highresdock
					delta_score = get_pose_ddg(working_fxn_, working_pose_);

					ms_tr.Debug << "Pre-move delta score = " << delta_score << ", fa_atr = " << fa_atr << ", fa_rep = " << fa_rep << ", fa_atr_rep before = " << fa_atr_rep_score_before << std::endl;

					//attempt to use movers to optimize placement a little more
					//begin setting up objects to make HighResDocker mover work

					//add ligand to pose residue type set if it is not already in the set
					//program doesn't work if the ligand isn't added to the pose residue type sets
					//also add it to the original pdb
					if ( ligand_added == false ) {
						ligand_added = true;

						//add the ligand mrt to the working and original poses
						add_ligand_to_pose_residuetypeset(lig_mrt);
					}

					//apply the highresdocker to working_pose
					run_HighResDock_on_working_pose(my_HighResDocker);

					//use working function to get ddg after highresdock
					delta_score = get_pose_ddg(working_fxn_, working_pose_);

					//this section up until collecting motifs off the placed ligand I do not think makes sense to put into its own function
					//My reasoning is that this could be made into 2 small specific scoring functions, which seems like extra work
					//The creation of the PDB comments uses a lot of upstream variables that would be a pain to pass into a function, so it may be better to just leave it inline in discover()

					//declaration of fa_atr_rep_score_after (may not be used)
					core::Real fa_atr_rep_score_after = 0;
					//optional check of fa_atr_rep after running highresdock
					//post_highresdock_fa_atr_rep_score
					if ( use_atr_rep ) {
						fa_atr_rep_score_after = fa_atr_rep_fxn_->score(*working_pose_);

						//check if after is worse than cutoff
						if ( fa_atr_rep_score_after > fa_atr_rep_cutoff_ ) {
							reset_working_pose();
							++clashing_counter;
							continue;
						}

						//reset atr and rep values based on score of whole pose
						fa_rep = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_rep];
						fa_atr = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_atr];
					}

					core::Real whole_score = 0;

					//run optional check of using ligand.wts score function as additional cutoff before attempting to keep
					if ( use_ligand_wts ) {
						// score the function with the whole ligand.tws function
						whole_score = whole_score_fxn_->score(*working_pose_);

						//check if the score is below the cutoff, kill if higher
						if ( whole_score > whole_fxn_cutoff_ ) {
							reset_working_pose();
							++clashing_counter;
							continue;
						}

						//print the score

						ms_tr.Debug << "ligand.wts score function score: " << whole_score << std::endl;

						//reset atr and rep values based on score of whole pose
						fa_rep = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_rep];
						fa_atr = working_pose_->energies().residue_total_energies(working_pose_->size())[core::scoring::fa_atr];
					}

					//only print post dock delta score if we did a score with the whole or atrrep function
					if ( use_atr_rep || use_ligand_wts ) {
						//ms_tr << "Post-dock delta score = " << delta_score << ", fa_atr = " << fa_atr << ", fa_rep = " << fa_rep << ", coordinate_constraint = " << sc_constraint_check << std::endl;

						ms_tr.Debug << "Post-dock delta score = " << delta_score << ", fa_atr = " << fa_atr << ", fa_rep = " << fa_rep;

						//with atr_rep
						if ( use_atr_rep ) {
							ms_tr.Debug << ", fa_atr_rep after = " << fa_atr_rep_score_after;
						}
						if ( use_ligand_wts ) {
							//without atr_rep
							ms_tr.Debug << ", whole = " << whole_score;
						}
						ms_tr.Debug << std::endl;
					}
					//check if whole_score is within cutoff, kill if not
					//need to remove ligand from poses so that they can be recycled
					//in theory, atr and rep should only improve, but this check helps make sure of that
					if ( delta_score > ddg_cutoff_ || fa_atr > fa_atr_cutoff_ || fa_rep > fa_rep_cutoff_ ) {
						reset_working_pose();
						++clashing_counter;
						continue;
					}

					//declare strings to serve as a small simplified table of information to be added to the placement pdb comments that can be easily extracted for later analysis
					//1 string will be a header line and the other will contain data that corresponds to the header
					//data will be in a csv format for easier parsing
					std::string comment_table_header = "";
					std::string comment_table_data = "";

					//name the pdb  that could come from the pose
					//current naming convention
					std::string pdb_name = output_prefix + "_ResPos_" + std::to_string(working_position) + "_ResID_" + discovery_position_residue + "_Trio" + std::to_string(i) + "_" + ligresOP->name() + "_motif_" + motifcop->remark();

					//add comments to working_pose for print
					//core::pose::add_comment(*working_pose_, "", );
					core::pose::add_comment(*working_pose_, "Placement: Output prefix:", output_prefix);
					core::pose::add_comment(*working_pose_, "Placement: Anchor residue index:", std::to_string(working_position));
					core::pose::add_comment(*working_pose_, "Placement: Anchor residue type:", discovery_position_residue);
					core::pose::add_comment(*working_pose_, "Placement: Ligand trio number:", std::to_string(i));
					core::pose::add_comment(*working_pose_, "Placement: Ligand name:", ligresOP->name());
					core::pose::add_comment(*working_pose_, "Placement: Placed motif remark:", motifcop->remark());

					//preceeding commas to account for the comment map keys
					comment_table_header = "ligand_name,ligand_atom_trio,source_pdb,anchor_residue_index,anchor_residue_type,placed_motif_remark,";
					comment_table_data = ligresOP->name() + "," + std::to_string(i) + "," + output_prefix + "," + std::to_string(working_position) + "," + discovery_position_residue + "," + motifcop->remark() + ",";

					//make a string that is the pdb name up to the motif that is used for motif collection of the placement (if that is used)
					std::string pdb_short_unique_name = pdb_name;

					//only use fa_rep and atr if we use atrrep or whole to score the whole system and pull how we fixed them
					if ( use_atr_rep || use_ligand_wts ) {
						pdb_name = pdb_name + "_rep_" + std::to_string(fa_rep) + "_atr_" + std::to_string(fa_atr);
						core::pose::add_comment(*working_pose_, "Scoring: Ligand fa_atr:", std::to_string(fa_atr));
						core::pose::add_comment(*working_pose_, "Scoring: Ligand fa_rep:", std::to_string(fa_rep));
						comment_table_data = comment_table_data + std::to_string(fa_atr) + "," + std::to_string(fa_rep) + ",";
					} else {
						//blank if we have no data, but keep placeholder to allow for easier mecshing of any data that did have the data
						comment_table_data = comment_table_data +  ","  + ",";
					}
					comment_table_header = comment_table_header + "fa_atr,fa_rep,ddg,";
					comment_table_data = comment_table_data + std::to_string(delta_score) + ",";

					//add delta (keeping in same order as I have in the past of atr and rep potentially being first)
					pdb_name = pdb_name + "_delta_" + std::to_string(delta_score);
					core::pose::add_comment(*working_pose_, "Scoring: Post-HighResDock system ddG:", std::to_string(delta_score));

					//adjust if using optional atr_rep post highresdock
					if ( use_atr_rep ) {
						pdb_name = pdb_name + "_atrrep_" + std::to_string(fa_atr_rep_score_after);
						core::pose::add_comment(*working_pose_, "Scoring: System fa_atr and fa_rep score:", std::to_string(fa_atr_rep_score_after));

						comment_table_data = comment_table_data + std::to_string(fa_atr_rep_score_after) + ",";
					} else {
						comment_table_data = comment_table_data + ",";
					}

					comment_table_header = comment_table_header + "fa_atr_rep,";

					//adjust file name if using ligand.wts
					if ( use_ligand_wts ) {
						pdb_name = pdb_name + "_whole_" + std::to_string(whole_score);
						core::pose::add_comment(*working_pose_, "Scoring: Whole system score:", std::to_string(whole_score));
						comment_table_data = comment_table_data + std::to_string(whole_score) + ",";
					} else {
						comment_table_data = comment_table_data + ",";
					}

					comment_table_header = comment_table_header + "whole_score,";

					//option to try to pull motifs from the passed placement and see what motifs are collected, how many there are, if motifs are made with any residues of interest, and if the motifs match any motifs in the motif library
					if ( option[ OptionKeys::motifs::collect_motifs_from_placed_ligand] ) {

						ms_tr.Trace << "Preparing to collect motif data from placed ligand." << std::endl;

						//make a vector that holds the following data in its indices as follows:
						//1 - total motifs made
						//2 - effectively a bool; 0 indicates that at least one pre-selected residue did not get a motif, 1 indicates that all did
						//3 - total motifs made against significant residues
						//4 - motifs that are considered close enough to at least one motif in the compare_library (considered real)
						//5 - ratio of motifs that are considered close enough compared to the total number of motifs that are collected (real motif ratio)
						utility::vector1< core::Real > placement_motifs_data;

						//run evaluate_motifs_of_pose from the ilm, which returns a vector with data on motifs from the placed ligand
						//comments are written to the pose, so we do not have to do that in this function now
						placement_motifs_data = ilm.evaluate_motifs_of_pose(working_pose_, mymap, pdb_short_unique_name);

						ms_tr.Trace << "Collected motif data from placed ligand." << std::endl;

						ms_tr.Debug << "Contents of placement_motifs_data vector: " << placement_motifs_data[1] << "," << placement_motifs_data[2] << "," << placement_motifs_data[3] << "," << placement_motifs_data[4] << "," << placement_motifs_data[5] << "," <<std::endl;

						//run checks with cutoffs
						//if minimum number of motifs made is not enough, kill placement
						//if at least one mandatory residue did not get a motif, kill
						//if the number of significant motifs made is greater than or equal to the cutoff, keep the placement, otherwise kill
						//determine if real ratio is greater than cutoff: minimum_ratio_of_real_motifs_from_ligand
						if ( placement_motifs_data[1] < min_motifs_cutoff_ || placement_motifs_data[2] == 0 || placement_motifs_data[3] < min_sig_motifs_cutoff_ || placement_motifs_data[5] < real_motif_ratio_cutoff_ ) {
							reset_working_pose();
							++clashing_counter;
							continue;
						}

						//work on building pdb name string if we passed and building the comment table strings
						pdb_name = pdb_name + "_motifs_" + std::to_string(placement_motifs_data[1]);
						comment_table_data = comment_table_data + std::to_string(placement_motifs_data[1]) + ",";

						if ( option[ OptionKeys::motifs::significant_residues_for_motifs].user() ) {
							pdb_name = pdb_name + "_sigmotifs_" + std::to_string(placement_motifs_data[3]);
							comment_table_data = comment_table_data + std::to_string(placement_motifs_data[3]) + ",";
						}
						comment_table_data = comment_table_data + ",";

						if ( option[ OptionKeys::motifs::check_if_ligand_motifs_match_real] ) {
							//add real motif ratio to pdb name
							pdb_name = pdb_name + "_realmotifratio_" + std::to_string(placement_motifs_data[5]);
							comment_table_data = comment_table_data + std::to_string(placement_motifs_data[4]) + "," + std::to_string(placement_motifs_data[5]) + ",";
						}
						comment_table_data = comment_table_data + ",,";

						//after optional modifications, add ".pdb" to cap off name
						pdb_name = pdb_name + ".pdb";

					} else {
						comment_table_data = comment_table_data + "," + "," + "," + ",";
					}

					comment_table_header = comment_table_header + "total_motif_count,significant_motif_count,real_motif_count,real_motif_ratio,";

					//after optional modifications, add ".pdb" to cap off name
					pdb_name = pdb_name + ".pdb";

					//if using clean naming scheme, replace name with clean name
					if ( option[ OptionKeys::motifs::clean_pdb_name] ) {
						pdb_name = output_prefix + "_" + ligresOP->name() + "_" + std::to_string(unique_placement_counter) + ".pdb";
						++unique_placement_counter;
					}

					//include the file name in the table to help with potential accession
					comment_table_data = "," + pdb_name + "," + comment_table_data;
					comment_table_header = ",pdb_file_name," + comment_table_header;

					//add comments for the table data
					core::pose::add_comment(*working_pose_, "Table Header", comment_table_header);
					core::pose::add_comment(*working_pose_, "Table Values", comment_table_data);

					//if we are keeping all placements that are better than a given ddg cutoff, it is better to not inflate the best_100_placements vector since we only use it for sorting
					//instead, just make the pdb and keep going
					if ( best_pdbs_to_keep == 0 ) {

						ms_tr.Debug << "Making pdb file for: " << pdb_name << std::endl;

						core::io::pdb::dump_pdb(*working_pose_, pdb_name);

						reset_working_pose();
						++passing_counter;

						//continue since we do not need to touch best_100_placements
						continue;
					}

					//if we are keeping a set number of placements, then we need to append to a list that will be sorted
					std::tuple<core::Real, core::pose::Pose, std::string> pose_tuple(delta_score, *working_pose_, pdb_name);

					reset_working_pose();

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


				ms_tr.Debug << "# passing cases for this trio = " << passing_counter << std::endl;
				ms_tr.Debug << "# clashing cases for this trio = " << clashing_counter << std::endl;
				ms_tr.Debug << "# post motif filter cases for this trio = " << pose_atom_check_counter << std::endl;


				ligand_passing_counter += passing_counter;
				ligand_clashing_counter += clashing_counter;


				ms_tr.Debug << "Total cases for trio: " << passing_counter + clashing_counter << std::endl;


			}//end of trio loop

			ms_tr << "Done iterating all trios, moving to next ligand" << std::endl;
			ms_tr << "Total passing attempts for ligand is " << ligand_passing_counter << std::endl;
			ms_tr << "Total clashing attempts for ligand is " << ligand_clashing_counter << std::endl;

			//assign the new cutoff once we hit 100 entries
			if ( best_placements.size() >= best_pdbs_to_keep && best_pdbs_to_keep != 0 ) {
				ddg_cutoff_ = std::get<0>(best_placements[0]);

				ms_tr.Debug << "New ddg cutoff is: " << ddg_cutoff_ << std::endl;

			}
		}
		//sort passing_placements if keeping only a specified amount
		if ( best_pdbs_to_keep != 0 ) {
			best_placements = mergeSort(best_placements);
		}

		//create pdbs of best scoring poses for pdb

		for ( const auto & best_pose_runner : best_placements ) {

			ms_tr.Debug << "Making pdb file for: " << std::get<2>(best_pose_runner) << std::endl;

			core::io::pdb::dump_pdb( ( std::get<1>(best_pose_runner)), std::get<2>(best_pose_runner));
		}


		ms_tr.Debug << "Number of placements that passed all filters: " << passed_placement_counter << std::endl;

	}
	return 1;
}

// @brief A function to reset the working_pose_ back to its original state, as found in original_pose
//to call this function the least amount of times (and minimize wasted time/resources), this function is only to be called right before a placement is culled before a continue statement after working_pose_ had a ligand added to it
void LigandDiscoverySearch::reset_working_pose()
{
	//assign working_pose_ to clone of original
	working_pose_ = original_pose_.clone();
}

//function to get a sub-library of motifs from the main library, based on the residue being used (only get for select residue)
//function may have use outside discover, so public can use it
protocols::motifs::MotifCOPs LigandDiscoverySearch::get_motif_sublibrary_by_aa(std::string residue_name)
{
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


	ms_tr.Debug << "Created motif sub-library for residue " << residue_name << " with " << motif_library_for_select_residue_.size() << " motifs in it." << std::endl;


	return motif_holder;
}
