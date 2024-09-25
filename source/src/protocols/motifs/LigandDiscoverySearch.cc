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
#include <map>
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

#include <protocols/motifs/IdentifyLigandMotifs.hh>

#include <utility/Binary_Util.hh>

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
	//make tracer for debugging
	static basic::Tracer ms_tr( "LigandDiscoverySearch.seed_cutoff_values", basic::t_info );

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

	//placement motif metric cutoffs
	ms_tr.Debug << "Using maximum placement motif to real motif RMSD distance: "  << option[ OptionKeys::motifs::duplicate_dist_cutoff] << std::endl;
	dist_threshold_ =  option[ basic::options::OptionKeys::motifs::duplicate_dist_cutoff ];

	ms_tr.Debug << "Using maximum placement motif to real motif angle distance: "  << option[ OptionKeys::motifs::duplicate_angle_cutoff] << std::endl;
	angl_threshold_ = option[ basic::options::OptionKeys::motifs::duplicate_angle_cutoff ];
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

	for (const auto & working_position : working_positions)
	{
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

// @brief function to set value of dist_threshold_
void LigandDiscoverySearch::set_distance_threshold(core::Size distance_threshold)
{
	dist_threshold_ = distance_threshold;
}

// @brief function to set value of angl_threshold_
void LigandDiscoverySearch::set_angl_threshold(core::Size angle_threshold)
{
	angl_threshold_ = angle_threshold;
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
		for ( const auto & atom_j :  atom_i_connects) {
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
	//make a tracer
	static basic::Tracer ms_tr( "LigandDiscoverySearch.make_minipose", basic::t_info );

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

	//check if fa_rep is good
	//positive score is bad
	//best scores are negative and closest to 0
	if ( fa_rep > fa_rep_cutoff_ ) {
		return false;
	}

	fa_atr_fxn_->score(*minipose);

	fa_atr = minipose->energies().residue_total_energies(minipose->size())[core::scoring::fa_atr];

	//run fa_atr check
	//don't keep if fa_atr is greater than cutoff
	if ( fa_atr > fa_atr_cutoff_ ) {
		return false;
	}

	//score whole minipose with atr_rep function
	fa_atr_rep_score_before = fa_atr_rep_fxn_->score(*minipose);

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
	sc_cst_set->add_constraint( core::scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_1, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_1 ), fx1 ) ) );

	core::scoring::func::FuncOP fx2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	sc_cst_set->add_constraint( core::scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_2, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_2 ), fx2 ) ) );

	core::scoring::func::FuncOP fx3( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	sc_cst_set->add_constraint( core::scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >( core::id::AtomID( trip_atom_3, working_pose_->size() ), core::id::AtomID( working_pose_->residue( working_position ).atom_index( "CA" ), 1 ), ligresOP->xyz( trip_atom_3 ), fx3 ) ) );

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

	//make a ResidueTypeSetCOP that will be pulled from the PoseResidueTypeSetOP
	core::chemical::ResidueTypeSetCOP def_rts(rts->default_rts());

	//use the name_mapOP function to get a residue type pointer based on the name of the ligand (could be either a nullptr or a pointer to the type)
	//should be a nullptr if it isn't in the set
	core::chemical::ResidueTypeCOP lig_rt(def_rts->name_mapOP(working_pose_->residue(working_pose_->size()).name()));

	//add the ligand to the residue type set if lig_rt is a nullpointer (which it should be)
	if ( lig_rt == nullptr ) {
		rts->add_base_residue_type(lig_mrt);
	}

	//reset the residue type sets for the working_pose_ and the original_pose_ so that they are updated to have the new ligand
	working_pose_->conformation().reset_residue_type_set_for_conf(rts);
	original_pose_.conformation().reset_residue_type_set_for_conf(rts);
}

//main function to run ligand discovery operations
//needs to have values set for working_pose_, motif_library_, and all_residues_
//parameter is a string to be a prefix name to use for outputted file names
core::Size LigandDiscoverySearch::discover(std::string output_prefix)
{
	//create tracer to identify points of the run
	//This tracer uses the following tracer outpurs: standard (no extension), Debug, Warning, and Trace
	static basic::Tracer ms_tr( "LigandDiscoverySearch.discover", basic::t_info );

	//create identifyligandmotifs object for use with deriving motifs from placed ligands
	IdentifyLigandMotifs ilm;

	// Make an atomtypeset to get atomtype integers for use in derive_adjacent_atoms_of_ligand
	core::chemical::AtomTypeSetCOP atset = core::chemical::ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );


	//iterate over all indices in working_positions_
	//if the size of working_positions is 0, return -1 because we want at least 1 index to work with
	if ( working_positions_.size() == 0 ) {

		ms_tr.Warning << "Length of working_positions_ vector is 0 and we have no positions to work with, killing the attempt now." << std::endl;

		return -1;
	}

	//run discovery over each position in working_positions_
	for (const auto & working_position : working_positions_) {

		ms_tr << "Current anchor residue position: " << working_position << std::endl;


		//declare empty discovery_position_residue name that will be properly written in initialization
		std::string discovery_position_residue = "";

		//determine whether to kill due to bad initialization
		bool kill_bad_init = false;

		//kill if working position is out of bounds
		if ( working_position < 1 || working_position > working_pose_->size() ) {

			ms_tr.Warning << "Working position of " << working_position << " is invalid as it is not a valid index to access a residue in your pose of size " << working_pose_->size() << std::endl;

			kill_bad_init = true;
		}
		else
		{
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
		std::map<std::tuple<std::string, std::string, std::string, std::string, std::string, std::string, std::string>,protocols::motifs::MotifCOPs> mymap;
		//only hash if we use the flag check_if_ligand_motifs_match_real, since this is currently the only thing we use it for
		if ( option[ OptionKeys::motifs::check_if_ligand_motifs_match_real] ) {
			hash_motif_library_into_map(motif_library_,mymap);
		}

		//set up target_residues_sf_ and target_residues_contact_
		//can interact with any of 3 integer_value vectors flags: in::target_residues, motifs::targer_residues_sf, and motifs::target_residues_contact
		//target_residues_sf will only write into the target_residues_sf vector; target_residues_contact will only write into target_residues_contact_
		//if target_residues flag is used, it will overwrite call to other flags and write the same contents to both vectors (easy way to fill both as same)
		//acceptable behavior (with handling) for either vector to be empty

		//read in target_residues_sf and target_residues_contact flags if they were used
		if ( option[ OptionKeys::motifs::target_residues_sf ].user() ) {
			target_residues_sf_ = option[ OptionKeys::motifs::target_residues_sf ]();
		}

		//read in target_residues and override
		if ( option[ OptionKeys::in::target_residues ].user() ) {
			target_residues_sf_ = option[ OptionKeys::in::target_residues ]();
		}

		core::Size x_shift = 0;
		core::Size y_shift = 0;
		core::Size z_shift = 0;
		int x_bound_int = 0;
		int y_bound_int = 0;
		int z_bound_int = 0;

		//create protein atom matrix
		create_protein_representation_matrix(x_shift, y_shift, z_shift, x_bound_int, y_bound_int, z_bound_int);

		//variable to  hold the linear scale by which we will increase the resolution of identification
		//default is 1. Probably want to set to 3-6
		//declaring here so it can be passed in calls for space fill function
		int resolution_increase_factor = option[ OptionKeys::motifs::resolution_scale_factor ];

		//declare variable set for space fill protein matrix
		//create vectors to hold values that correspond to xyz values in the space fill matrix
		//these probably shouldn't be xyzVector objects, since they don't directly correspond to xyz coordinates, but rather that they are values that correspond to xyz axes, scalars, and segments

		//used to note the shift of the vector so that the smallest xyz values correspond to 0,0,0
		utility::vector1<core::Size>  xyz_shift_sf (3,0);

		//used to note the maxima of the xyz vector
		utility::vector1<core::Size>  xyz_bound_sf (3,0);

		//used to note the minimum values of points in the sub-area of the space fill vector
		utility::vector1<core::Size>  sub_xyz_min_sf (3,0);

		//used to note the maximum values of points in the sub-area of the space fill vector
		utility::vector1<core::Size>  sub_xyz_max_sf (3,0);

		//create 2 vectors to hold occupied ratio scores (core::Real) and raw counts of occupied, unoccupied, and total cells (core::Size) for both the main area and sub area

		//index 1 corresponds to the occupied ratio for the whole system
		//index 2 corresponds to only the sub-area
		utility::vector1<core::Real>  occupied_ratios (2,0);

		//whole system
		//index 1 = occupied cell count
		//index 2 = unoccupied cell count
		//index 3 = total cell count
		//sub area
		//index 4 = occupied cell count
		//index 5 = unoccpuied cell count
		//index 6 = total cell count
		utility::vector1<core::Size>  matrix_data_counts (6,0);

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

		//create the space fill representation matrix
		//do even if no cutoff scores were made so that the information can be generated and returned anyway
		create_protein_representation_matrix_space_fill(xyz_shift_sf, xyz_bound_sf, resolution_increase_factor, sub_xyz_min_sf, sub_xyz_max_sf, occupied_ratios, matrix_data_counts, working_position);

		//create a dummy mutableresiduetype that is empty (all atoms deleted) so we can use it for the expore_space_fill_matrix function
		core::conformation::ResidueOP dummyligresOP = all_residues_[1];

		//derive a mutable residue type
		core::chemical::MutableResidueTypeOP dummylig_mrt( new core::chemical::MutableResidueType( dummyligresOP->type() ) );

		//run through each atom in dummylig_mrt and delete it so it is effectively clear
		//this definitely goes against what the data structure is designed for (considering that the default constructor is blocked from use), but hopefully this just does what I need it to do
		//see if this works: keep deleting the atom at index 1 until the number of atoms in the mrt is 0
		while ( dummylig_mrt->natoms() > 0 )
				{
			//debug print of the first atom in the list of verticies and number of atoms

			ms_tr.Trace << "Number of atoms in dummy mrt is now: " <<  dummylig_mrt->natoms() << std::endl;

			//pop the first atom in the list of vertices
			dummylig_mrt->delete_atom(dummylig_mrt->all_atoms()[1]);
		}


		ms_tr.Trace << "Number of atoms in dummy mrt is now: " <<  dummylig_mrt->natoms() << std::endl;

		//delete all chis in the dummy ligand
		//probably easiest to use delete_terminal_chi
		while ( dummylig_mrt->nchi() > 0 )
				{

					ms_tr.Trace << "Number of chi in residue is currently: " << dummylig_mrt->nchi() << std::endl;

			//delete the terminal chi
			dummylig_mrt->delete_terminal_chi();
		}

		ms_tr.Trace << "Number of chi in residue is currently: " << dummylig_mrt->nchi() << std::endl;


		//create a pdb of the space fill matrix if the flag selects for it
		if ( option[ OptionKeys::motifs::output_space_fill_matrix_pdbs ] ) {
			//call function to export, use "empty" as the string prefix
			export_space_fill_matrix_as_C_H_O_N_pdb(protein_representation_matrix_space_fill_, xyz_shift_sf, xyz_bound_sf, resolution_increase_factor,
				occupied_ratios, "empty", *dummylig_mrt);
			//export_space_fill_matrix_as_C_H_O_N_pdb(protein_representation_matrix_, xyz_shift_sf, xyz_bound_sf, resolution_increase_factor,
			// sub_xyz_min_sf, sub_xyz_max_sf, occupied_ratios, "empty");
		}

		//make a copy of matrix_data_counts that represents the system without a placed ligand
		//this will be used so that we can modify the original as we look at systems with ligands placed, and can revert to the original
		utility::vector1<core::Size>  matrix_data_counts_empty = matrix_data_counts;

		//for whole ligand.wts function, determine also if we want to use the function for scoring at all
		bool use_ligand_wts = option[ OptionKeys::motifs::score_with_ligand_wts_function ];

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

		//create a copy of the working_pose_ that working_pose_ can be reset to after each placement attempt
		original_pose_(*working_pose_);

		//hold the number of placements that pass all filters and could enter the top 100 placements
		int passed_placement_counter = 0;

		//now we have the filtered motif library to work with, run through each  atom trio in each ligand and try to match it against all motifs for the residue
		for (const auto & tracker : all_residues_ ) {

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
					has_clashing = ligand_clash_check(ligresOP, x_shift, y_shift, z_shift, x_bound_int, y_bound_int, z_bound_int);

					//continue because we clash
					if ( has_clashing == true ) {
						continue;
					}

					//space fill analysis block
					//faster than using score functions
					//this is used to determine if enough of a defined binding pocked is filled by a placed ligand or not (ensures elimination of off-target placements)
					//only attempt this block if we request at least one of the 2 score cutoffs for space filling

					if ( option[ OptionKeys::motifs::space_fill_cutoff_score ].user() || option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() || option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {

						//track the space fill score for the system and sub_area, hold as a 2d Real vector
						//index 1 is the score for the whole system, index 2 is the score for the sub area
						//default values are 0
						utility::vector1<core::Real> space_fill_scores (2,0);

						//run space fill analysis function, set to a vector
						SpaceFillMatrix current_space_fill_matrix = space_fill_analysis(ligresOP, xyz_shift_sf, xyz_bound_sf, resolution_increase_factor, sub_xyz_min_sf, sub_xyz_max_sf, space_fill_scores, matrix_data_counts);

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
							export_space_fill_matrix_as_C_H_O_N_pdb(current_space_fill_matrix, xyz_shift_sf, xyz_bound_sf, resolution_increase_factor,
								space_fill_scores, matrix_pdb_prefix, *dummylig_mrt);
						}

						//at end before check, reset matrix_data_counts so that it returns to the empty state
						matrix_data_counts = matrix_data_counts_empty;

						//run check for if the placement is passable based on score
						//either the whole system score or the sub area score need to pass (one can fail)
						if ( (space_fill_scores[1] >= score_cutoffs_sf[1] and option[ OptionKeys::motifs::space_fill_cutoff_score ].user()) || (space_fill_scores[2] >= score_cutoffs_sf[2] and option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user()) || (fill_differential >= score_cutoffs_sf[3] and option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user()) ) {
							//print out where any scores passed

							if ( space_fill_scores[1] >= score_cutoffs_sf[1] and option[ OptionKeys::motifs::space_fill_cutoff_score ].user() ) {
								ms_tr.Debug << space_fill_scores[1] << " : whole system, passed" << std::endl;
							} else if ( option[ OptionKeys::motifs::space_fill_cutoff_score ].user() ) {
								ms_tr.Debug << space_fill_scores[1] << " : whole system, failed" << std::endl;
							}

							if ( space_fill_scores[2] >= score_cutoffs_sf[2] and option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() ) {
								ms_tr.Debug << space_fill_scores[2] << " : sub-area, passed" << std::endl;
							} else if ( option[ OptionKeys::motifs::space_fill_cutoff_score_sub ].user() ) {
								ms_tr.Debug << space_fill_scores[2] << " : sub-area, failed" << std::endl;
							}

							if ( fill_differential >= score_cutoffs_sf[3] and option[ OptionKeys::motifs::space_fill_cutoff_differential_score_sub ].user() ) {
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
						if (make_minipose(minipose, ligresOP) == false)
						{
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
					if ( !minipose_scoring)
					{
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

					//use fa atr/rep or whole function based on highresdock_with_whole_score_fxn flag to get ddg before highresdock
					if ( option[ OptionKeys::motifs::highresdock_with_whole_score_fxn ] ) {
						//whole
						delta_score = get_pose_ddg(whole_score_fxn_, working_pose_);
					} else {
						//atrrep
						delta_score = get_pose_ddg(fa_atr_rep_fxn_, working_pose_);
					}
					
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
					protocols::ligand_docking::HighResDockerOP my_HighResDocker_atrrep( new protocols::ligand_docking::HighResDocker(3, 0, fa_atr_rep_fxn_, my_movemapbuilder) );
					protocols::ligand_docking::HighResDockerOP my_HighResDocker_whole( new protocols::ligand_docking::HighResDocker(3, 0, whole_score_fxn_, my_movemapbuilder) );


					//set highresdockers to not repack
					my_HighResDocker_atrrep->set_allow_repacking(false);
					my_HighResDocker_whole->set_allow_repacking(false);

					//apply with highresdocker for desired score function
					if ( option[ OptionKeys::motifs::highresdock_with_whole_score_fxn ] ) {
						my_HighResDocker_atrrep->apply(*working_pose_);
					} else {
						my_HighResDocker_whole->apply(*working_pose_);
					}

					//use fa atr/rep or whole function based on highresdock_with_whole_score_fxn flag to get new ddg
					if ( option[ OptionKeys::motifs::highresdock_with_whole_score_fxn ] ) {
						//whole
						delta_score = get_pose_ddg(whole_score_fxn_, working_pose_);
					} else {
						//atrrep
						delta_score = get_pose_ddg(fa_atr_rep_fxn_, working_pose_);
					}

					//declaration of fa_atr_rep_score_after (may not be used)
					core::Real fa_atr_rep_score_after = 0;
					//optional check of fa_atr_rep after running highresdock
					//post_highresdock_fa_atr_rep_score
					if ( option[ OptionKeys::motifs::post_highresdock_fa_atr_rep_score ] ) {
						//whole_score_fxn_->score(*working_pose_);
						//5/1/24 replacing score with the atr_rep function instead of whole
						fa_atr_rep_score_after = fa_atr_rep_fxn_->score(*working_pose_);

						//check if after is worse than cutoff
						if ( fa_atr_rep_score_after > fa_atr_rep_cutoff_ ) {
							working_pose_->delete_residue_slow(working_pose_->size());
							//create new poseop of the original pose (to wipe any highresdock changes) to the pose
							core::pose::PoseOP original_poseop(original_pose_.clone());
							//assign the original_poseop to working_pose_
							working_pose_ = original_poseop;
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
							//reset the working pose because its whole score was not good enough
							working_pose_->delete_residue_slow(working_pose_->size());
							//create new poseop of the original pose (to wipe any highresdock changes) to the pose
							core::pose::PoseOP original_poseop(original_pose_.clone());
							//assign the original_poseop to working_pose_
							working_pose_ = original_poseop;
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
					if (option[ OptionKeys::motifs::post_highresdock_fa_atr_rep_score ] || use_ligand_wts ) {
						//ms_tr << "Post-dock delta score = " << delta_score << ", fa_atr = " << fa_atr << ", fa_rep = " << fa_rep << ", coordinate_constraint = " << sc_constraint_check << std::endl;

						ms_tr.Debug << "Post-dock delta score = " << delta_score << ", fa_atr = " << fa_atr << ", fa_rep = " << fa_rep;

						//with atr_rep
						if ( option[ OptionKeys::motifs::post_highresdock_fa_atr_rep_score ] ) {
							ms_tr.Debug << ", fa_atr_rep after = " << fa_atr_rep_score_after;
						} else if ( use_ligand_wts ) {
							//without atr_rep
							ms_tr.Debug << ", whole = " << whole_score;
						}
						ms_tr.Debug << std::endl;
					}
					//check if whole_score is within cutoff, kill if not
					//need to remove ligand from poses so that they can be recycled
					//in theory, atr and rep should only improve, but this check helps make sure of that
					if ( delta_score > ddg_cutoff_ || fa_atr > fa_atr_cutoff_ || fa_rep > fa_rep_cutoff_ ) {
						working_pose_->delete_residue_slow(working_pose_->size());
						//create new poseop of the original pose (to wipe any highresdock changes) to the pose
						core::pose::PoseOP original_poseop(original_pose_.clone());
						//assign the original_poseop to working_pose_
						working_pose_ = original_poseop;
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
					if ( option[ OptionKeys::motifs::post_highresdock_fa_atr_rep_score ] || use_ligand_wts ) {
						pdb_name = pdb_name + "_rep_" + std::to_string(fa_rep) + "_atr_" + std::to_string(fa_atr);
						core::pose::add_comment(*working_pose_, "Scoring: Ligand fa_atr:", std::to_string(fa_atr));
						core::pose::add_comment(*working_pose_, "Scoring: Ligand fa_rep:", std::to_string(fa_rep));
						comment_table_data = comment_table_data + std::to_string(fa_atr) + "," + std::to_string(fa_rep) + ",";
					}
					else
					{
						//blank if we have no data, but keep placeholder to allow for easier mecshing of any data that did have the data
						comment_table_data = comment_table_data +  ","  + ",";
					}
					comment_table_header = comment_table_header + "fa_atr,fa_rep,ddg,";
					comment_table_data = comment_table_data + std::to_string(delta_score) + ",";

					//add delta (keeping in same order as I have in the past of atr and rep potentially being first)
					pdb_name = pdb_name + "_delta_" + std::to_string(delta_score);
					core::pose::add_comment(*working_pose_, "Scoring: Post-HighResDock system ddG:", std::to_string(delta_score));

					//adjust if using optional atr_rep post highresdock
					if ( option[ OptionKeys::motifs::post_highresdock_fa_atr_rep_score ] ) {
						pdb_name = pdb_name + "_atrrep_" + std::to_string(fa_atr_rep_score_after);
						core::pose::add_comment(*working_pose_, "Scoring: System fa_atr and fa_rep score:", std::to_string(fa_atr_rep_score_after));

						comment_table_data = comment_table_data + std::to_string(fa_atr_rep_score_after) + ",";
					}
					else
					{
						comment_table_data = comment_table_data + ",";
					}

					comment_table_header = comment_table_header + "fa_atr_rep,";

					//adjust file name if using ligand.wts
					if ( use_ligand_wts ) {
						pdb_name = pdb_name + "_whole_" + std::to_string(whole_score);
						core::pose::add_comment(*working_pose_, "Scoring: Whole system score:", std::to_string(whole_score));
						comment_table_data = comment_table_data + std::to_string(whole_score) + ",";
					}
					else
					{
						comment_table_data = comment_table_data + ",";
					}

					comment_table_header = comment_table_header + "whole_score,";

					//option to try to pull motifs from the passed placement and see what motifs are collected, how many there are, if motifs are made with any residues of interest, and if the motifs match any motifs in the motif library
					if ( option[ OptionKeys::motifs::collect_motifs_from_placed_ligand] ) {
						//create a new motif library to hold motifs
						protocols::motifs::MotifLibrary placement_library;

						//make vector that holds the indices of residues that contribute to motifs (probably the easiest way to track if motifs were made on residues of interest)
						utility::vector1< Size > prot_pos_that_made_motifs_size;

						//use ilm process_for_motifs to obtain motifs from the pose
						ilm.process_for_motifs(*working_pose_, pdb_short_unique_name, placement_library, prot_pos_that_made_motifs_size);

						//convery the prot_pos vector from size to int (easier to use int because this interacts with vectors that are pulled from args that don't seem to be able to be pulled as size; I can convert those to size, but that is extra work)
						utility::vector1< int > prot_pos_that_made_motifs = prot_pos_that_made_motifs_size;

						//determine how many motifs were made and how many were made on significant residues
						core::Size motifs_made = prot_pos_that_made_motifs.size();

						ms_tr.Debug << "Ligand placement created " << motifs_made << " total motifs" << std::endl;						

						//if minimum number of motifs made is not enough, kill placement
						if ( motifs_made < min_motifs_cutoff_ ) {
							working_pose_->delete_residue_slow(working_pose_->size());
							//create new poseop of the original pose (to wipe any highresdock changes) to the pose
							core::pose::PoseOP original_poseop(original_pose_.clone());
							//assign the original_poseop to working_pose_
							working_pose_ = original_poseop;
							++clashing_counter;
							continue;
						}

						pdb_name = pdb_name + "_motifs_" + std::to_string(motifs_made);
						core::pose::add_comment(*working_pose_, "Placement motifs: Total motifs made:", std::to_string(motifs_made));

						comment_table_data = comment_table_data + std::to_string(motifs_made) + ",";

						//convert the motif library to motifCOPS
						protocols::motifs::MotifCOPs placement_libraryCOPs = placement_library.library();

						//iterate over the library of motifs created by the ligand placement to add a comment for each that is the motif's remark
						//create a counter too to print out each unique motif
						core::Size placement_motif_counter = 0;
						for ( auto ligmotifcop : placement_libraryCOPs ) {

							core::pose::add_comment(*working_pose_, "Placement motifs: Placement motif " + std::to_string(placement_motif_counter) + ":", ligmotifcop->remark());

							++placement_motif_counter;
						}

						//check if there are motifs made for all mandatory residues
						if ( option[ OptionKeys::motifs::mandatory_residues_for_motifs].user() ) {
							//bool to help control loops to determine whether to kill the placed ligand
							bool kill = false;
							utility::vector1< int > mandatory_residues_for_motifs = option[ OptionKeys::motifs::mandatory_residues_for_motifs] ;
							for ( core::Size sig_res_pos = 1; sig_res_pos < mandatory_residues_for_motifs.size(); ++sig_res_pos ) {
								//kill unless we get a match of a motif made having the same residue index as the current residue in the mandatory list
								kill = true;
								for ( core::Size motif_made = 1; motif_made < prot_pos_that_made_motifs.size(); ++motif_made ) {
									if ( prot_pos_that_made_motifs[motif_made] == mandatory_residues_for_motifs[sig_res_pos] ) {
										//tick up the counter for significant motifs made if there is a match in the residue index for the motif and a significant residue
										kill = false;
									}
								}

								//if kill is still true, we didn't get a motif for the mandatory residue, move forward with killing the ligand
								if ( kill ) {
									
									ms_tr.Debug << "Not motifs made for residue index " << sig_res_pos << std::endl;
									
									break;
								}
							}
							if ( kill ) {
								working_pose_->delete_residue_slow(working_pose_->size());
								//create new poseop of the original pose (to wipe any highresdock changes) to the pose
								core::pose::PoseOP original_poseop(original_pose_.clone());
								//assign the original_poseop to working_pose_
								working_pose_ = original_poseop;
								++clashing_counter;
								continue;
							}

							
							ms_tr.Trace << "Made motifs for all mandatory residues" << std::endl;
							

						}

						//variable to track how many motifs hit a significant residue
						core::Size significant_motifs_made = 0;

						if ( option[ OptionKeys::motifs::significant_residues_for_motifs].user() ) {
							
							ms_tr.Trace << "Ligand placement created motifs against significant residues: ";
							

							//build a string that holds the significant indices that had a motif form against it
							std::string significant_residue_string = "";

							utility::vector1< int > significant_residues_for_motifs = option[ OptionKeys::motifs::significant_residues_for_motifs] ;
							for ( core::Size sig_res_pos = 1; sig_res_pos < significant_residues_for_motifs.size(); ++sig_res_pos ) {
								for ( core::Size motif_made = 1; motif_made < prot_pos_that_made_motifs.size(); ++motif_made ) {
									if ( prot_pos_that_made_motifs[motif_made] == significant_residues_for_motifs[sig_res_pos] ) {
										//tick up the counter for significant motifs made if there is a match in the residue index for the motif and a significant residue
										++significant_motifs_made;

										ms_tr.Debug << prot_pos_that_made_motifs[motif_made] << ",";
										significant_residue_string += std::to_string(prot_pos_that_made_motifs[motif_made]) + ",";
									}
								}
							}

							
							ms_tr.Debug << std::endl;
							ms_tr.Debug << "Ligand placement created " << significant_motifs_made << " motifs for significant residues" << std::endl;
							

							pdb_name = pdb_name + "_sigmotifs_" + std::to_string(significant_motifs_made);

							core::pose::add_comment(*working_pose_, "Placement motifs: Motifs made against significant residues count:", std::to_string(significant_motifs_made));
							core::pose::add_comment(*working_pose_, "Placement motifs: Motifs made against significant residues:", significant_residue_string);

							comment_table_data = comment_table_data + std::to_string(significant_motifs_made) + ",";
						}
						else
						{
							comment_table_data = comment_table_data + ",";
						}

						//if the number of significant motifs made is greater than or equal to the cutoff, keep the placement, otherwise kill
						if ( significant_motifs_made < min_sig_motifs_cutoff_ ) {
							working_pose_->delete_residue_slow(working_pose_->size());
							//create new poseop of the original pose (to wipe any highresdock changes) to the pose
							core::pose::PoseOP original_poseop(original_pose_.clone());
							//assign the original_poseop to working_pose_
							working_pose_ = original_poseop;
							++clashing_counter;
							continue;
						}

						//check if motifs that were generated match real motifs (as inputted into this program as the motif list)
						//also determine if the ratio of real generated motifs is above the expected cutoff
						if ( option[ OptionKeys::motifs::check_if_ligand_motifs_match_real] ) {
							//counter to count the number of motifs generated by the placement that are close enough to a real motif we have
							//use to derive a ratio (with potential for placement to be filtered with a cutoff)
							core::Size motifs_that_look_real = 0;

							//iterate over the library of motifs created by the ligand placement
							for ( auto ligmotifcop : placement_libraryCOPs ) {
								//derive tuple key
								std::tuple<std::string, std::string, std::string, std::string, std::string, std::string, std::string> curkey_tuple(ligmotifcop->restype_name1(),ligmotifcop->res1_atom1_name(),ligmotifcop->res1_atom2_name(),ligmotifcop->res1_atom3_name(),ligmotifcop->res2_atom1_name(),ligmotifcop->res2_atom2_name(),ligmotifcop->res2_atom3_name());

								//pull out the motif library that matches the current motif that we are on by atom names (if there is one)
								//use map count function to determine if the key exists
								if ( mymap.count(curkey_tuple) > 0 ) {
									//key exists
									//pull out motif library at key address and then compare all motifs in the list against ligmotifcop to see if it resembles a real motif
									protocols::motifs::MotifCOPs real_motifs = mymap[curkey_tuple];

									//create a bool that determines if we found a real match for the motif or not
									bool real_match_found = false;

									//iterate over the library and compare
									for ( auto realmotifcop : real_motifs ) {
										//compare code based on remove_duplicate_motifs
										//difference is that we are not deleting any motifs, instead if we get a match hit, we will call out ligand-derived motif real
										//secondary confirm that the motifs match (there is no good reason why we should hit a continue off of this if the map is formed right)
										if ( ligmotifcop->restype_name1() != realmotifcop->restype_name1() ) continue;
										if ( ligmotifcop->res1_atom1_name() != realmotifcop->res1_atom1_name() ) continue;
										if ( ligmotifcop->res1_atom2_name() != realmotifcop->res1_atom2_name() ) continue;
										if ( ligmotifcop->res1_atom3_name() != realmotifcop->res1_atom3_name() ) continue;
										if ( ligmotifcop->res2_atom1_name() != realmotifcop->res2_atom1_name() ) continue;
										if ( ligmotifcop->res2_atom2_name() != realmotifcop->res2_atom2_name() ) continue;
										if ( ligmotifcop->res2_atom3_name() != realmotifcop->res2_atom3_name() ) continue;

										core::Real motif_distance = 0;
										core::Real motif_theta = 0;

										jump_distance(ligmotifcop->forward_jump(), realmotifcop->forward_jump(), motif_distance, motif_theta);
										
										ms_tr.Debug << "Comparing to motif " << realmotifcop->remark() << std::endl;
										ms_tr.Debug << "Distance: " << motif_distance << std::endl;
										ms_tr.Debug << "Theta: " << motif_theta << std::endl;
										

										if ( motif_distance < dist_threshold_ && motif_theta < angl_threshold_ ) {
											//note that the motif matches a real one
											
											ms_tr.Debug << "Current motif matches real motif with distance: " << motif_distance << "  and angle difference: "  << motif_theta << std::endl;
											if ( ligmotifcop->has_remark() ) {
												ms_tr.Debug << "Ligand remark is: " << ligmotifcop->remark() << std::endl;
											}
											if ( realmotifcop->has_remark() ) {
												ms_tr.Debug << "Real motif remark is: " << realmotifcop->remark() << std::endl;
											}
											

											//increment real counter
											++motifs_that_look_real;

											std::string real_motif_match_info = "remark: " + realmotifcop->remark() + ", distance: " + std::to_string(motif_distance) + ", angle: " + std::to_string(motif_theta); 

											//add comment about real motif match for placement motif
											core::pose::add_comment(*working_pose_, "Placement motifs: Real motif check - " + ligmotifcop->remark(), real_motif_match_info);

											real_match_found = true;

											//greedy algorithm, stop for current ligmotifcop when we hit the first match since we got what we wanted
											break;
										}

									}

									//if no real match was found, note comment
									if (real_match_found == false){
									core::pose::add_comment(*working_pose_, "Placement motifs: Real motif check - " + ligmotifcop->remark(), "No real match");
									}
								} else {
									//key (and real motif in our library) does not exist
									//declare that this motif is not considered real
									
									ms_tr.Debug << "No real motifs identified for motif: " << *ligmotifcop << std::endl;									

									core::pose::add_comment(*working_pose_, "Placement motifs: Real motif check - " + ligmotifcop->remark(), "No real match");
								}
							}

							//determine the ratio of ligand motifs that match real ones
							//convert numerator an denominator to real to ensure we have no floating point issues
							core::Real motifs_that_look_real_real = motifs_that_look_real;
							core::Real motifs_made_real = motifs_made;
							core::Real real_looking_motif_ratio = motifs_that_look_real_real/motifs_made_real;

							//determine if real ratio is greater than cutoff: minimum_ratio_of_real_motifs_from_ligand
							if ( real_looking_motif_ratio < real_motif_ratio_cutoff_ ) {
								//kill placement
								working_pose_->delete_residue_slow(working_pose_->size());
								//create new poseop of the original pose (to wipe any highresdock changes) to the pose
								core::pose::PoseOP original_poseop(original_pose_.clone());
								//assign the original_poseop to working_pose_
								working_pose_ = original_poseop;
								++clashing_counter;
								continue;
							}

							//add real motif ratio to pdb name
							pdb_name = pdb_name + "_realmotifratio_" + std::to_string(real_looking_motif_ratio);

							core::pose::add_comment(*working_pose_, "Placement motifs: Real motif count:", std::to_string(motifs_that_look_real_real));
							core::pose::add_comment(*working_pose_, "Placement motifs: Real motif ratio:", std::to_string(real_looking_motif_ratio));
							comment_table_data = comment_table_data + std::to_string(motifs_that_look_real_real) + "," + std::to_string(real_looking_motif_ratio) + ",";

						}
						else
						{
							comment_table_data = comment_table_data + "," + ",";
						}

					}
					else
					{
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

						//directly remove residue from end of pose
						working_pose_->delete_residue_slow(working_pose_->size());
						//create new poseop of the original pose (to wipe any highresdock changes) to the pose
						core::pose::PoseOP original_poseop(original_pose_.clone());
						//assign the original_poseop to working_pose_
						working_pose_ = original_poseop;
						++passing_counter;

						//continue since we do not need to touch best_100_placements
						continue;
					}

					//if we are keeping a set number of placements, then we need to append to a list that will be sorted
					std::tuple<core::Real, core::pose::Pose, std::string> pose_tuple(delta_score, *working_pose_, pdb_name);

					//directly remove residue from end of pose
					working_pose_->delete_residue_slow(working_pose_->size());
					//create new poseop of the original pose (to wipe any highresdock changes) to the pose
					core::pose::PoseOP original_poseop(original_pose_.clone());
					//assign the original_poseop to working_pose_
					working_pose_ = original_poseop;

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

//function to get a sub-library of motifs from the main library, based on the residue being used (only get for select residue)
//function may have use outside discover, so public can use it
protocols::motifs::MotifCOPs LigandDiscoverySearch::get_motif_sublibrary_by_aa(std::string residue_name)
{
	static basic::Tracer ms_tr( "LigandDiscoverySearch.get_motif_sublibrary_by_aa", basic::t_info );

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

//create protein_representation_matrix_
//uses working_pose to make the matrix
void LigandDiscoverySearch::create_protein_representation_matrix(core::Size & x_shift, core::Size & y_shift, core::Size & z_shift, int & x_bound_int, int & y_bound_int, int & z_bound_int)
{

	//create tracer to identify points of the run
	static basic::Tracer ms_tr( "LigandDiscoverySearch.create_protein_representation_matrix", basic::t_info );

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
	
	ms_tr.Debug << "Creating protein clash coordinate matrix. Dimensions of matrix are " << x_bound << "," << y_bound << "," << z_bound << std::endl;
	

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

//create protein_representation_matrix_space_fill_
//uses working_pose to make the matrix
//void LigandDiscoverySearch::create_protein_representation_matrix_space_fill(core::Size & x_shift, core::Size & y_shift, core::Size & z_shift, core::Size & x_bound_int, core::Size & y_bound_int, core::Size & z_bound_int, int & resolution_increase_factor,
//core::Size & sub_x_min, core::Size & sub_x_max, core::Size & sub_y_min, core::Size & sub_y_max, core::Size & sub_z_min, core::Size & sub_z_max, core::Real & occupied_ratio, core::Real & sub_occupied_ratio)
//condensing arguments in function to use vectors to hold xyz trios
void LigandDiscoverySearch::create_protein_representation_matrix_space_fill(utility::vector1<core::Size> & xyz_shift, utility::vector1<core::Size> & xyz_bound, int & resolution_increase_factor,
	utility::vector1<core::Size> & sub_xyz_min, utility::vector1<core::Size> & sub_xyz_max, utility::vector1<core::Real> & occupied_ratios, utility::vector1<core::Size> & matrix_data_counts, core::Size working_position)
{
	//values that can be seeded into the matrix with their meaning. Only 0-3 can be seeded using this function, and then space_fill_analysis can seed 4+
	//note, 0 and 2 are only even numbers that are currently used for now
	/*
	0 = empty and out of sub area, carbon, black
	1 = protein and out of sub area, fluorine, icy blue
	2 = empty and in sub area, oxygen, red
	3 = protein and in sub area, nitrogen, blue
	4 = do not use, keep even for unoccupied space
	5 = ligand and out of sub area, sulphur, yellow
	6 = do not use, keep even for unoccupied space
	7 = ligand and in sub area, chlorine, green
	8 = do not use, keep even for unoccupied space
	9 = ligand and protein out of sub area, phosphorous, orange
	10 = do not use, keep even for unoccupied space
	11 = ligand and protein in sub area, iodine, purple
	*/

	//create tracer to identify points of the run
	static basic::Tracer ms_tr( "LigandDiscoverySearch.create_protein_matrix_space_fill", basic::t_info );

	int smallest_x = 1;
	int smallest_y = 1;
	int smallest_z = 1;

	int largest_x = 1;
	int largest_y = 1;
	int largest_z = 1;

	//create a list of coordinates of each atom to hold and work with to fill the protein_representation_matrix
	//can't seem to make a vector of xyzVector objects, so will need to just make a custome 2D vector  to  hold the data
	//unlike clash detection, going to use floats as well for increaserd precision
	utility::vector1<numeric::xyzVector<int>> atom_coordinates;

	//vector to hold the coordinates of atoms and the lennard jones radius
	utility::vector1<utility::vector1<core::Real>> atom_coordinates_float_and_lj_radius;

	//determine largest and smallest x,y,z  values to determine dimensions of matrix
	for ( core::Size res_num = 1; res_num <= working_pose_->size(); ++res_num ) {
		for ( core::Size atom_num = 1; atom_num <= working_pose_->residue(res_num).natoms(); ++atom_num ) {
			//get the x,y,z data of the atom
			numeric::xyzVector<int> atom_xyz = working_pose_->residue(res_num).xyz(atom_num);

			utility::vector1<core::Real> atom_xyz_float_with_lj_radius;
			atom_xyz_float_with_lj_radius.push_back(atom_xyz.x());
			atom_xyz_float_with_lj_radius.push_back(atom_xyz.y());
			atom_xyz_float_with_lj_radius.push_back(atom_xyz.z());
			atom_xyz_float_with_lj_radius.push_back(working_pose_->residue(res_num).atom_type(atom_num).lj_radius());

			atom_coordinates_float_and_lj_radius.push_back(atom_xyz_float_with_lj_radius);


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

			//determine if any are the largest
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
	xyz_shift[1] = (smallest_x * -1) + 1;
	xyz_shift[2] = (smallest_y * -1) + 1;
	xyz_shift[3] = (smallest_z * -1) + 1;

	//apply shift values to largest to get boundaries
	xyz_bound[1] = xyz_shift[1] + largest_x;
	xyz_bound[2] = xyz_shift[2] + largest_y;
	xyz_bound[3] = xyz_shift[3] + largest_z;

	//apply constant shifts to all coordinates
	//apply resolution factor to all atoms
	for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates.size(); ++xyzVec ) {
		atom_coordinates[xyzVec].x() += xyz_shift[1];
		atom_coordinates[xyzVec].y() += xyz_shift[2];
		atom_coordinates[xyzVec].z() += xyz_shift[3];

		atom_coordinates[xyzVec].x() *= resolution_increase_factor;
		atom_coordinates[xyzVec].y() *= resolution_increase_factor;
		atom_coordinates[xyzVec].z() *= resolution_increase_factor;

		atom_coordinates_float_and_lj_radius[xyzVec][1] += xyz_shift[1];
		atom_coordinates_float_and_lj_radius[xyzVec][2] += xyz_shift[2];
		atom_coordinates_float_and_lj_radius[xyzVec][3] += xyz_shift[3];

		atom_coordinates_float_and_lj_radius[xyzVec][1] *= resolution_increase_factor;
		atom_coordinates_float_and_lj_radius[xyzVec][2] *= resolution_increase_factor;
		atom_coordinates_float_and_lj_radius[xyzVec][3] *= resolution_increase_factor;
		atom_coordinates_float_and_lj_radius[xyzVec][4] *= resolution_increase_factor;
	}

	//create 3D matrix to roughly represent 3D coordinate space of protein
	//bad syntax, may just have to do  an iterative  fill

	xyz_bound[1] *= resolution_increase_factor;
	xyz_bound[2] *= resolution_increase_factor;
	xyz_bound[3] *= resolution_increase_factor;

	
	ms_tr.Debug << "Creating space fill matrix. Dimensions of matrix are " << xyz_bound[1] << "," << xyz_bound[2] << "," << xyz_bound[3] << std::endl;


	//make a matrix, and we will copy it over to the global one once we make it
	SpaceFillMatrix protein_representation_matrix;

	for ( core::Size x = 1; x <= xyz_bound[1]; ++x ) {
		//make a 2D  matrix
		utility::vector1<utility::vector1<core::Size>> sub_matrix;

		for ( core::Size y = 1; y <= xyz_bound[2]; ++y ) {

			//make a 1D matrix, seed with 0 values
			utility::vector1<core::Size> sub_sub_matrix(xyz_bound[3],  0);
			//push 1D  matrix into 2D
			sub_matrix.push_back(sub_sub_matrix);

		}
		//push a 2D  matrix into the 3D matrix
		protein_representation_matrix.push_back(sub_matrix);
	}

	//seed the matrix with approximate coordinates of each atom
	//use the vector that has the LJ radii to fill out cells
	for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates_float_and_lj_radius.size(); ++xyzVec ) {
		//iterate through all cells in the matrix and determine if the sphere projected by the atom center and LJ radius hits this cell
		//there may be a more efficient way to do this...
		//I think there is! instead of iterating the whole volume of the matrix, only iterate about a range bound by a cube with with side length = 2 x LJ radius
		//adjust min and max to ensure we don't run off of the matrix

		core::Size x_min = atom_coordinates_float_and_lj_radius[xyzVec][1] - atom_coordinates_float_and_lj_radius[xyzVec][4];
		core::Size x_max = atom_coordinates_float_and_lj_radius[xyzVec][1] + atom_coordinates_float_and_lj_radius[xyzVec][4];
		if ( x_min < 1 ) {
			x_min = 1;
		}
		if ( x_max > xyz_bound[1] ) {
			x_max = xyz_bound[1];
		}
		core::Size y_min = atom_coordinates_float_and_lj_radius[xyzVec][2] - atom_coordinates_float_and_lj_radius[xyzVec][4];
		core::Size y_max = atom_coordinates_float_and_lj_radius[xyzVec][2] + atom_coordinates_float_and_lj_radius[xyzVec][4];
		if ( y_min < 1 ) {
			y_min = 1;
		}
		if ( y_max > xyz_bound[2] ) {
			y_max = xyz_bound[2];
		}
		core::Size z_min = atom_coordinates_float_and_lj_radius[xyzVec][3] - atom_coordinates_float_and_lj_radius[xyzVec][4];
		core::Size z_max = atom_coordinates_float_and_lj_radius[xyzVec][3] + atom_coordinates_float_and_lj_radius[xyzVec][4];
		if ( z_min < 1 ) {
			z_min = 1;
		}
		if ( z_max > xyz_bound[3] ) {
			z_max = xyz_bound[3];
		}

		//iterate over each cell between xmin-max, ymin-max, zmin-max, a cube around the atom
		for ( core::Size x = x_min; x <= x_max; ++x ) {
			for ( core::Size y = y_min; y <= y_max; ++y ) {
				for ( core::Size z = z_min; z <= z_max; ++z ) {
					//use distance formula to figure out if cell x,y,z is within the sphere projected by the atom point about its LJ radius
					//get distance between x,y,z and the atom point
					core::Real atom_cell_distance = sqrt(((x - atom_coordinates_float_and_lj_radius[xyzVec][1]) * (x - atom_coordinates_float_and_lj_radius[xyzVec][1])) + ((y - atom_coordinates_float_and_lj_radius[xyzVec][2]) * (y - atom_coordinates_float_and_lj_radius[xyzVec][2])) + ((z - atom_coordinates_float_and_lj_radius[xyzVec][3]) * (z - atom_coordinates_float_and_lj_radius[xyzVec][3])));
					//if distance is less than the radius, then the point is occupied
					if ( atom_cell_distance < atom_coordinates_float_and_lj_radius[xyzVec][4] ) {
						protein_representation_matrix[x][y][z] = 1;
					}
				}
			}
		}
	}

	//system should be processed. Get a count of occupied vs unoccupied cells
	core::Real total_cells = xyz_bound[1] * xyz_bound[2] * xyz_bound[3];

	core::Real occupied_cell_count = 0;
	core::Real unoccupied_cell_count = 0;

	for ( core::Size x = 1; x <= xyz_bound[1]; ++x ) {
		for ( core::Size y = 1; y <= xyz_bound[2]; ++y ) {
			for ( core::Size z = 1; z <= xyz_bound[3]; ++z ) {
				if ( protein_representation_matrix[x][y][z] == 1 ) {
					++occupied_cell_count;
				} else {
					++unoccupied_cell_count;
				}
			}
		}
	}

	occupied_ratios[1] = occupied_cell_count / total_cells;

	
	ms_tr.Debug << "Total: " << total_cells << std::endl;
	ms_tr.Debug << "Occupied: " << occupied_cell_count << std::endl;
	ms_tr.Debug << "Unoccupied: " << unoccupied_cell_count << std::endl;
	ms_tr.Debug << "Occupied-Total Ratio: " << occupied_ratios[1] << std::endl;
	

	//write values to matrix_data_counts
	matrix_data_counts[1] = occupied_cell_count;
	matrix_data_counts[2] = unoccupied_cell_count;
	matrix_data_counts[3] = total_cells;

	//NEW FOR BINDING POCKET FILL TEST
	//Use the anchor residue and inputted vector of residues of interest (if listed), and create a box to use to represent a binding pocket (subsection of whole pose) for small scale space fill analysis

	//define subsection if no vector of residues of interest listed (based on NBR of anchor * 2, similar to minipose derivation)
	//get location of anchor residue center/NBR atom
	//#########

	//check to see if binding_pocket_center_sf flag is used; create a vector to correspond to its coordinates later when defining the max and min for the sub-area
	bool using_binding_pocket_center_sf = false;
	if ( option[ OptionKeys::motifs::binding_pocket_center_sf ].user() ) {
		using_binding_pocket_center_sf = true;
	}

	core::Real sub_total_cells = 0;
	core::Real sub_occupied_cell_count = 0;
	core::Real sub_unoccupied_cell_count = 0;

	//define the sub-area min and max for xyz
	sub_xyz_min[1] = 0;
	sub_xyz_max[1] = 0;
	sub_xyz_min[2] = 0;
	sub_xyz_max[2] = 0;
	sub_xyz_min[3] = 0;
	sub_xyz_max[3] = 0;

	//default of only looking at area about anchor residue; probably not the best method if you actually care about using this method
	if ( target_residues_sf_.size() == 0 && using_binding_pocket_center_sf == false ) {

		
		ms_tr.Trace << "Defining sub-area only by anchor residue." << std::endl;
		

		//define nbr_atom_xyz vector as the xyz vector of the nbr atom of the xyz residue
		numeric::xyzVector<int> nbr_atom_xyz = working_pose_->residue(working_position).nbr_atom_xyz();

		//get the nbr radius of the residue
		core::Real anchor_nbr_radius = working_pose_->residue(working_position).nbr_radius();

		//debugging: get coordinates of nbr atom and nbr radius
		
		ms_tr.Debug << "Raw coordinates of nbr atom and nbr radius: " << nbr_atom_xyz.x() << ", " << nbr_atom_xyz.y() << ", " << nbr_atom_xyz.z() << "; " << anchor_nbr_radius << std::endl;
		

		//apply the xyz shift
		nbr_atom_xyz[0] += static_cast<int>(xyz_shift[1]);
		nbr_atom_xyz[1] += static_cast<int>(xyz_shift[2]);
		nbr_atom_xyz[2] += static_cast<int>(xyz_shift[3]);

		
		ms_tr.Debug << "Shifted not scaled coordinates of nbr atom and nbr radius: " << nbr_atom_xyz[0] << ", " << nbr_atom_xyz[1] << ", " << nbr_atom_xyz[2] << "; " << anchor_nbr_radius << std::endl;
		
		//translate NBR radius length to define the sub-area to investigate
		nbr_atom_xyz[0] *= resolution_increase_factor;
		nbr_atom_xyz[1] *= resolution_increase_factor;
		nbr_atom_xyz[2] *= resolution_increase_factor;

		anchor_nbr_radius *= resolution_increase_factor;

		//define the sub-area min and max for xyz
		sub_xyz_min[1] = nbr_atom_xyz[0] - anchor_nbr_radius;
		sub_xyz_max[1] = nbr_atom_xyz[0] + anchor_nbr_radius;
		sub_xyz_min[2] = nbr_atom_xyz[1] - anchor_nbr_radius;
		sub_xyz_max[2] = nbr_atom_xyz[1] + anchor_nbr_radius;
		sub_xyz_min[3] = nbr_atom_xyz[2] - anchor_nbr_radius;
		sub_xyz_max[3] = nbr_atom_xyz[2] + anchor_nbr_radius;


		//debugging: get coordinates of nbr atom and nbr radius
		
		ms_tr.Debug << "xyz shifts and resolution factor: " << xyz_shift[1] << ", " << xyz_shift[2] << ", " << xyz_shift[3] <<  "; " << resolution_increase_factor << std::endl;
		ms_tr.Debug << "Shifted and scaled coordinates of nbr atom and nbr radius: " << nbr_atom_xyz[0] << ", " << nbr_atom_xyz[1] << ", " << nbr_atom_xyz[2] << "; " << anchor_nbr_radius << std::endl;
		ms_tr.Debug << "Min and Max xyz: " << sub_xyz_min[1] << ", " << sub_xyz_min[2] << ", " << sub_xyz_min[3] << "; " << sub_xyz_max[1] << ", " << sub_xyz_max[2] << ", " << sub_xyz_max[3] << " " << std::endl;
		

	} else if ( target_residues_sf_.size() > 0 ) {
		//
		//define subsection if the vector of residues is listed
		//if target_residues_sf and binding_pocket_center_sf flags are used, the cube/sphere area will be defined by the min and max values from both. This can help make sure that specific residue regions are covered (could turn the cube into a rectangular prism and sphere into an elipsoid)
		ms_tr.Trace << "Defining sub-area by inputted residues using motifs::target_residues_sf flag." << std::endl;
		

		//run through the vector of indices and identify the coordinates of the nbr atom
		//could technically check all atoms in each residue, but would be longer and nbr should get us a good enough estimate
		for ( const auto & target_residue : target_residues_sf_ ) {
			//get the scaled xyz coordinates of the residue's nbr atom and adjust the sub xyz minmax values
			//have as int to account for negative coordinates (which we will convert to positive by shifting)
			numeric::xyzVector<int> nbr_atom_xyz_int = working_pose_->residue(target_residue).nbr_atom_xyz();

			//apply the xyz shift to each value (and do before applying the resolution increase factor)
			nbr_atom_xyz_int[0] += xyz_shift[1];
			nbr_atom_xyz_int[1] += xyz_shift[2];
			nbr_atom_xyz_int[2] += xyz_shift[3];

			//translate NBR radius length to define the sub-area to investigate
			nbr_atom_xyz_int[0] *= resolution_increase_factor;
			nbr_atom_xyz_int[1] *= resolution_increase_factor;
			nbr_atom_xyz_int[2] *= resolution_increase_factor;

			//create new xyzvector of the nbr atom xyz that is shifted and scaled
			numeric::xyzVector<core::Size> nbr_atom_xyz = nbr_atom_xyz_int;

			// get min and max for x,y,z for subarea

			//x
			if ( nbr_atom_xyz[0] < sub_xyz_min[1] || sub_xyz_min[1] == 0 ) {

				sub_xyz_min[1] = nbr_atom_xyz[0];

			}
			if ( nbr_atom_xyz[0] > sub_xyz_max[1] || sub_xyz_max[1] == 0 ) {
				//we check later if we exceed the max, so not a worry here
				sub_xyz_max[1] = nbr_atom_xyz[0];
			}

			//y
			if ( nbr_atom_xyz[1] < sub_xyz_min[2] || sub_xyz_min[2] == 0 ) {

				sub_xyz_min[2] = nbr_atom_xyz[1];

			}
			if ( nbr_atom_xyz[1] > sub_xyz_max[2] || sub_xyz_max[2] == 0 ) {
				//we check later if we exceed the max, so not a worry here
				sub_xyz_max[2] = nbr_atom_xyz[1];
			}

			//z
			if ( nbr_atom_xyz[2] < sub_xyz_min[3] || sub_xyz_min[3] == 0 ) {

				sub_xyz_min[3] = nbr_atom_xyz[2];

			}
			if ( nbr_atom_xyz[2] > sub_xyz_max[3] || sub_xyz_max[3] == 0 ) {
				//we check later if we exceed the max, so not a worry here
				sub_xyz_max[3] = nbr_atom_xyz[2];
			}


		}
	} else if ( using_binding_pocket_center_sf == true ) {
		//pull the xyz coordinates of the binding pocket center
		//may the vector to hold the xyz coordinates
		utility::vector1<core::Size> binding_pocket_center_xyz = option[ OptionKeys::motifs::binding_pocket_center_sf ]();

		//make sure that the user has inputted 3 valid coordinates and flag a warning for each value beyond the first 3
		//center coordinates are invalid if they are outside the space fill matrix (apply the xyz shift and resolution scaling and then check)
		//iterate over each coordinate to apply the respectful xyz shift and resolution and make these checks as well as for bad coordinates
		//also ensure that there are at least 3 values. If there are not, we can not use the 0-2 coordinate values and we will kill the attempt at this method

		bool vector_big_enough = true;

		if ( binding_pocket_center_xyz.size() <= 2 ) {
			vector_big_enough = false;
		}

		if ( vector_big_enough ) {
			for ( core::Size coordinate = 1; coordinate <= binding_pocket_center_xyz.size(); ++coordinate ) {
				//apply relevant xyz shift as long as the coordinate is between 1-3, otherwise throw a warning
				if ( coordinate <= 3 ) {
					binding_pocket_center_xyz[coordinate] += xyz_shift[coordinate];
					binding_pocket_center_xyz[coordinate] *= resolution_increase_factor;
				} else {
					
					ms_tr.Warning << "You have inputted an excess value using the motifs::binding_pocket_center_sf flag. This is value #" << coordinate << " in the input vector with a value of: " << binding_pocket_center_xyz[coordinate] << ". Ignoring this bad input, and you should review the proper usage for this flag." << std::endl;
					
				}
			}
		} else {
			
			ms_tr.Warning << "You have only inputted " << binding_pocket_center_xyz.size() << " coordinates for the binding pocket center xyz coordinates. We need exactly 3 to work. We can't do anything with this, and will skip this method and will just use an area around the anchor residue." << std::endl;
			
		}

		//get the sub min and max values
		//need to test and make sure that potential underflow with Size data type doesn't mess things up
		//get the radius value
		//this is for a single value to be used for xyz
		core::Real radius_real = option[ OptionKeys::motifs::binding_pocket_radius_sf ]();
		//scale the radius by the resolution factor
		radius_real *= resolution_increase_factor;

		//declare a vector to use for x,y,z
		//index 0 corresponds to x, index 1 corresponds to y, index 2 corresponds to z
		utility::vector1<core::Real> binding_pocket_dimensions;
		binding_pocket_dimensions.push_back(radius_real);
		binding_pocket_dimensions.push_back(radius_real);
		binding_pocket_dimensions.push_back(radius_real);

		ms_tr << "binding_pocket_dimensions: " << binding_pocket_dimensions[1] << "," << binding_pocket_dimensions[2] << "," << binding_pocket_dimensions[3] << std::endl;
		

		//overwrite if we used binding_pocket_dimensions_sf flag (which can use a different value in x,y,z)
		if ( option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].user() ) {
			//make sure dimensions vector is valid
			//can not use if fewer than 3 values (do not use), warn if there are more than 3
			if ( option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size()<3 ) {
			ms_tr.Warning << "Input vector for binding_pocket_dimensions_sf is too small with only " << option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size() << " entries. We will default to using binding_pocket_radius_sf." << std::endl;
				
			} else {
				if ( option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size()>3 ) {
					
					ms_tr.Warning << "Input vector for binding_pocket_dimensions_sf has more than 3 entries with having " << option[ OptionKeys::motifs::binding_pocket_dimensions_sf ].size() << " entries. We will only use the first 3 values." << std::endl;
					
				}

				//set the first 3 values of the binding_pocket_dimensions_sf vector to binding_pocket_dimensions (and also apply the resolution factor)
				binding_pocket_dimensions = option[ OptionKeys::motifs::binding_pocket_dimensions_sf ]();

				//apply resolution factor to each dimension
				for ( core::Size dimension = 1; dimension <= binding_pocket_dimensions.size(); ++dimension ) {
					binding_pocket_dimensions[dimension] *= resolution_increase_factor;
				}

			}

		}


		ms_tr << "binding_pocket_dimensions: " << binding_pocket_dimensions[1] << "," << binding_pocket_dimensions[2] << "," << binding_pocket_dimensions[3] << std::endl;
		

		//apply the radius
		//make sure that we don't underflow for any values, and set the value to 1 if we would
		//values could underflow if the coordinate is really far in the negative direction
		//the shape will be junk if the entire shape is out of bounds
		//x
		if ( binding_pocket_center_xyz[1] - binding_pocket_dimensions[1] <= 0 ) {
			sub_xyz_min[1] = 1;
		} else {
			sub_xyz_min[1] = binding_pocket_center_xyz[1] - binding_pocket_dimensions[1];
		}
		if ( binding_pocket_center_xyz[1] + binding_pocket_dimensions[1] <= 0 ) {
			sub_xyz_max[1] = 1;
		} else {
			sub_xyz_max[1] = binding_pocket_center_xyz[1] + binding_pocket_dimensions[1];
		}

		//y
		if ( binding_pocket_center_xyz[2] - binding_pocket_dimensions[2] <= 0 ) {
			sub_xyz_min[2] = 1;
		} else {
			sub_xyz_min[2] = binding_pocket_center_xyz[2] - binding_pocket_dimensions[2];
		}
		if ( binding_pocket_center_xyz[2] + binding_pocket_dimensions[2] <= 0 ) {
			sub_xyz_max[2] = 1;
		} else {
			sub_xyz_max[2] = binding_pocket_center_xyz[2] + binding_pocket_dimensions[2];
		}

		//z
		if ( binding_pocket_center_xyz[3] - binding_pocket_dimensions[3] <= 0 ) {
			sub_xyz_min[3] = 1;
		} else {
			sub_xyz_min[3] = binding_pocket_center_xyz[3] - binding_pocket_dimensions[3];
		}
		if ( binding_pocket_center_xyz[3] + binding_pocket_dimensions[3] <= 0 ) {
			sub_xyz_max[3] = 1;
		} else {
			sub_xyz_max[3] = binding_pocket_center_xyz[3] + binding_pocket_dimensions[3];
		}
	}

	//if any values are outside the bounds of the main vector, adjust them to match the boundary and avoid going out of bounds
	if ( sub_xyz_min[1] < 1 ) {
		sub_xyz_min[1] = 1;
	}
	if ( sub_xyz_max[1] > xyz_bound[1] ) {
		sub_xyz_max[1] = xyz_bound[1];
	}

	if ( sub_xyz_min[2] < 1 ) {
		sub_xyz_min[2] = 1;
	}
	if ( sub_xyz_max[2] > xyz_bound[2] ) {
		sub_xyz_max[2] = xyz_bound[2];
	}

	if ( sub_xyz_min[3] < 1 ) {
		sub_xyz_min[3] = 1;
	}
	if ( sub_xyz_max[3] > xyz_bound[3] ) {
		sub_xyz_max[3] = xyz_bound[3];
	}

	
	ms_tr << "Sub-region stats without placed region, bound by adjusted coordinates: x(" << sub_xyz_min[1] << "->" << sub_xyz_max[1] << ") y(" << sub_xyz_min[2] << "->" << sub_xyz_max[2] << ") z(" << sub_xyz_min[3] << "->" << sub_xyz_max[3] << ")" << std::endl;
	
	//run through the boundaries and determine the space fill ratio of the sub-area
	for ( core::Size x = sub_xyz_min[1]; x <= sub_xyz_max[1]; ++x ) {
		for ( core::Size y = sub_xyz_min[2]; y <= sub_xyz_max[2]; ++y ) {
			for ( core::Size z = sub_xyz_min[3]; z <= sub_xyz_max[3]; ++z ) {
				if ( protein_representation_matrix[x][y][z] == 1 ) {
					++sub_occupied_cell_count;

					//adjust the value to now be 3 to indicate that this is full within the sub-area
					protein_representation_matrix[x][y][z] = 3;
				} else {
					++sub_unoccupied_cell_count;

					//adjust the value to now be 2 to indicate that this is empty within the sub-area
					protein_representation_matrix[x][y][z] = 2;
				}
				++sub_total_cells;
			}
		}
	}

	//calculate the ratio (to be passed by reference out of the function)
	occupied_ratios[2] = sub_occupied_cell_count / sub_total_cells;

	//output information on sub-region

	
	ms_tr.Debug << "Total: " << sub_total_cells << std::endl;
	ms_tr.Debug << "Occupied: " << sub_occupied_cell_count << std::endl;
	ms_tr.Debug << "Unoccupied: " << sub_unoccupied_cell_count << std::endl;
	ms_tr.Debug << "Occupied-Total Ratio: " << occupied_ratios[2] << std::endl;
	

	matrix_data_counts[4] = sub_occupied_cell_count;
	matrix_data_counts[5] = sub_unoccupied_cell_count;
	matrix_data_counts[6] = sub_total_cells;

	//copy the generated matrix over to the global variable
	protein_representation_matrix_space_fill_ = protein_representation_matrix;
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

// function to check if the space-filling capacity of the placed ligand is adequate
//returns a boolean value based on whether or not the space filling is satisfactory
//uses the protein representation space fill matrix
//
//satisfaction is based on whether the ratio of occupied cells to unoccupied cells is >= a user-inputted threshold
//returns space_fill_matrix_copy (the filled matrix with the ligand)

LigandDiscoverySearch::SpaceFillMatrix LigandDiscoverySearch::space_fill_analysis(core::conformation::ResidueOP ligresOP, utility::vector1<core::Size> & xyz_shift, utility::vector1<core::Size> & xyz_bound, int & resolution_increase_factor,
	utility::vector1<core::Size> & sub_xyz_min, utility::vector1<core::Size> & sub_xyz_max, utility::vector1<core::Real> & occupied_ratios, utility::vector1<core::Size> & matrix_data_counts)
{
	//values that can be seeded into the matrix with their meaning. Only 5,7,9,11 can be seeded using this function since we adjust values based on the ligand being present
	//note, 0 and 2 are only even numbers that are currently used for now
	/*
	0 = empty and out of sub area, carbon, black
	1 = protein and out of sub area, fluorine, icy blue
	2 = empty and in sub area, oxygen, red
	3 = protein and in sub area, nitrogen, blue
	4 = do not use, keep even for unoccupied space
	5 = ligand and out of sub area, sulphur, yellow
	6 = do not use, keep even for unoccupied space
	7 = ligand and in sub area, chlorine, green
	8 = do not use, keep even for unoccupied space
	9 = ligand and protein out of sub area, phosphorous, orange
	10 = do not use, keep even for unoccupied space
	11 = ligand and protein in sub area, iodine, purple
	*/

	//create tracer to identify points of the run
	static basic::Tracer ms_tr( "LigandDiscoverySearch.space_fill_analysis", basic::t_info );

	//debugging
	//print out what matrix data counts looks like before and after modification
	
	ms_tr.Debug << "matrix_data_counts: " << matrix_data_counts[1] << ", " << matrix_data_counts[2] << ", " << matrix_data_counts[3] << ", " << matrix_data_counts[4] << ", " << matrix_data_counts[5] << ", " << matrix_data_counts[6] << std::endl;

	//iterate through all atoms in the placed ligand, and update the count of occupied cells in the space fill matrix, going atom by atom
	//We will not edit the matrix, and will instead just update the count if a cells occupied by an atom in the ligand is not already occupied by an atom of the protein system

	//make a temporary copy of protein_representation_matrix_space_fill_ so that we can adjust the boolean values (and make sure that we can properly mark occupied cells once as we check with the atoms)
	//otherwise, we risk overcounting occupied cells (which also leads to underflow in the unoccupied cells)
	//make sure it is not a deep copy
	SpaceFillMatrix space_fill_matrix_copy = protein_representation_matrix_space_fill_;


	core::Size num_atoms_in_ligand = ligresOP->natoms();

	for ( core::Size residue_atom_iterator = 1; residue_atom_iterator <= num_atoms_in_ligand; ++residue_atom_iterator ) {
		//convert coordinates of current atom into format that can be read into space fill matrix matrix
		numeric::xyzVector<int> test_atom_xyz_int = ligresOP->xyz(residue_atom_iterator);

		//apply the xyz shift to each value (and do before applying the resolution increase factor)
		test_atom_xyz_int[0] += xyz_shift[1];
		test_atom_xyz_int[1] += xyz_shift[2];
		test_atom_xyz_int[2] += xyz_shift[3];

		//translate NBR radius length to define the sub-area to investigate
		test_atom_xyz_int[0] *= resolution_increase_factor;
		test_atom_xyz_int[1] *= resolution_increase_factor;
		test_atom_xyz_int[2] *= resolution_increase_factor;

		//temporary for debugging; list the coordinates of the atom before and after shifting
		//ms_tr << "Atom #" << residue_atom_iterator << " raw xyz: " << test_atom_xyz_int[0] << "," << test_atom_xyz_int[1] << "," << test_atom_xyz_int[2] << " adjusted xyz: " << test_atom_xyz_int[0] << "," << test_atom_xyz_int[1] << "," << test_atom_xyz_int[2] << std::endl;

		//create new xyzvector of the nbr atom xyz that is shifted and scaled
		numeric::xyzVector<core::Size> test_atom_xyz = test_atom_xyz_int;

		//get the lj radius of the atom and scale it with the resolution factor
		//atom_xyz_float_with_lj_radius.push_back(working_pose_->residue(res_num).atom_type(atom_num).lj_radius());
		//not storing in a vector this time, since we don't need it after this iteration of the for loop
		core::Real test_atom_lj_radius = ligresOP->atom_type(residue_atom_iterator).lj_radius() * resolution_increase_factor;

		//define the minimum and maximum areas to explore, based on the coordinate of the atom and the lj radius
		//ensure we don't run off the area

		core::Size x_min = test_atom_xyz_int[0] - test_atom_lj_radius;
		core::Size x_max = test_atom_xyz_int[0] + test_atom_lj_radius;
		if ( x_min < 1 ) {
			x_min = 1;
		}
		if ( x_max > xyz_bound[1] ) {
			x_max = xyz_bound[1];
		}
		core::Size y_min = test_atom_xyz_int[1] - test_atom_lj_radius;
		core::Size y_max = test_atom_xyz_int[1] + test_atom_lj_radius;
		if ( y_min < 1 ) {
			y_min = 1;
		}
		if ( y_max > xyz_bound[2] ) {
			y_max = xyz_bound[2];
		}
		core::Size z_min = test_atom_xyz_int[2] - test_atom_lj_radius;
		core::Size z_max = test_atom_xyz_int[2] + test_atom_lj_radius;
		if ( z_min < 1 ) {
			z_min = 1;
		}
		if ( z_max > xyz_bound[3] ) {
			z_max = xyz_bound[3];
		}

		//iterate over the xyz of this area
		//iterate over each cell between xmin-max, ymin-max, zmin-max, a cube around the atom
		for ( core::Size x = x_min; x <= x_max; ++x ) {
			for ( core::Size y = y_min; y <= y_max; ++y ) {
				for ( core::Size z = z_min; z <= z_max; ++z ) {

					//use distance formula to figure out if cell x,y,z is within the sphere projected by the atom point about its LJ radius
					//get distance between x,y,z and the atom point
					//distance = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
					core::Real atom_cell_distance = sqrt(((x - test_atom_xyz_int[0]) * (x - test_atom_xyz_int[0])) + ((y - test_atom_xyz_int[1]) * (y - test_atom_xyz_int[1])) + ((z - test_atom_xyz_int[2]) * (z - test_atom_xyz_int[2])));

					//if distance is less than the radius, then the cell/point is occupied by this atom and we should investigate
					if ( atom_cell_distance < test_atom_lj_radius ) {
						//make sure that the cell isn't occupied already in the space fill matrix
						//if not already occupied, adjust the occupied and unoccupied counts
						if ( space_fill_matrix_copy[x][y][z] == 0 ) {
							//std::cout << x << "," << y << "," << z << std::endl;
							//occupied by ligand and outside of sub area
							space_fill_matrix_copy[x][y][z] = 5;

							//adjust matrix data counts for occupied and unoccupied for the main system
							++matrix_data_counts[1];
							--matrix_data_counts[2];

							//if x,y,z lies within the boundary of the sub area, adjust the matrix data counts for it too
							if ( x >= sub_xyz_min[1] && x <= sub_xyz_max[1] && y >= sub_xyz_min[2] && y <= sub_xyz_max[2] && z >= sub_xyz_min[3] && z <= sub_xyz_max[3] ) {
								//adjust matrix data counts for occupied and unoccupied for the sub area
								++matrix_data_counts[4];
								--matrix_data_counts[5];
							}
						} else if ( space_fill_matrix_copy[x][y][z] == 1 ) {
							//std::cout << x << "," << y << "," << z << std::endl;
							//occupied by ligand and protein out of sub area
							space_fill_matrix_copy[x][y][z] = 9;

							//adjust matrix data counts for occupied and unoccupied for the main system
							++matrix_data_counts[1];
							--matrix_data_counts[2];

							//if x,y,z lies within the boundary of the sub area, adjust the matrix data counts for it too
							if ( x >= sub_xyz_min[1] && x <= sub_xyz_max[1] && y >= sub_xyz_min[2] && y <= sub_xyz_max[2] && z >= sub_xyz_min[3] && z <= sub_xyz_max[3] ) {
								//adjust matrix data counts for occupied and unoccupied for the sub area
								++matrix_data_counts[4];
								--matrix_data_counts[5];
							}
						} else if ( space_fill_matrix_copy[x][y][z] == 2 ) {
							//std::cout << x << "," << y << "," << z << std::endl;
							//occupied by ligand in sub area
							space_fill_matrix_copy[x][y][z] = 7;

							//adjust matrix data counts for occupied and unoccupied for the main system
							++matrix_data_counts[1];
							--matrix_data_counts[2];

							//if x,y,z lies within the boundary of the sub area, adjust the matrix data counts for it too
							if ( x >= sub_xyz_min[1] && x <= sub_xyz_max[1] && y >= sub_xyz_min[2] && y <= sub_xyz_max[2] && z >= sub_xyz_min[3] && z <= sub_xyz_max[3] ) {
								//adjust matrix data counts for occupied and unoccupied for the sub area
								++matrix_data_counts[4];
								--matrix_data_counts[5];
							}
						} else if ( space_fill_matrix_copy[x][y][z] == 3 ) {
							//std::cout << x << "," << y << "," << z << std::endl;
							//occupied by ligand and protein in sub area
							space_fill_matrix_copy[x][y][z] = 11;

							//adjust matrix data counts for occupied and unoccupied for the main system
							++matrix_data_counts[1];
							--matrix_data_counts[2];

							//if x,y,z lies within the boundary of the sub area, adjust the matrix data counts for it too
							if ( x >= sub_xyz_min[1] && x <= sub_xyz_max[1] && y >= sub_xyz_min[2] && y <= sub_xyz_max[2] && z >= sub_xyz_min[3] && z <= sub_xyz_max[3] ) {
								//adjust matrix data counts for occupied and unoccupied for the sub area
								++matrix_data_counts[4];
								--matrix_data_counts[5];
							}
						}
					}
				}
			}
		}
	}

	
	ms_tr.Debug << "matrix_data_counts: " << matrix_data_counts[1] << ", " << matrix_data_counts[2] << ", " << matrix_data_counts[3] << ", " << matrix_data_counts[4] << ", " << matrix_data_counts[5] << ", " << matrix_data_counts[6] << std::endl;
	
	//derive space fill scores for whole and sub areas
	occupied_ratios[1] = static_cast<core::Real>(matrix_data_counts[1])/static_cast<core::Real>(matrix_data_counts[3]);
	occupied_ratios[2] = static_cast<core::Real>(matrix_data_counts[4])/static_cast<core::Real>(matrix_data_counts[6]);
	return space_fill_matrix_copy;
}


// @brief debugging function to export a space fill matrix as a pdb. Occupied cells are represented as a nitrogen and unoccupied cells are represented as an oxygen (considering making one for a clash matrix too)
// if printing the whole matrix and not just the sub-area, occupied cells are represented by hydrogens and unoccupied are represented by carbon
//if desired, returns the pose created
//highly recommended to only use this for debugging and small-scale purposes, as this method is extremely slow compared to the rest of the protocol
pose::Pose LigandDiscoverySearch::export_space_fill_matrix_as_C_H_O_N_pdb(SpaceFillMatrix space_fill_matrix, utility::vector1<core::Size> & xyz_shift, utility::vector1<core::Size> & xyz_bound, int & resolution_increase_factor,
	utility::vector1<core::Real> & occupied_ratios, std::string pdb_name_prefix, core::chemical::MutableResidueType dummylig_mrt)
{
	//create tracer to identify points of the run
	static basic::Tracer ms_tr( "LigandDiscoverySearch.export_space_fill_matrix_as_C_H_O_N_pdb", basic::t_info );

	std::string matrix_pdb_name = pdb_name_prefix + "_WholeRatio_" + std::to_string(occupied_ratios[1]) + "_SubRatio_" + std::to_string(occupied_ratios[2]) + ".pdb";
	
	ms_tr.Trace << "Preparing to make visualization pose for " << matrix_pdb_name << std::endl;
	

	//create a new mutableresiduetype to work on adding atoms to
	//base it from the dummy mrt that had all of its atoms deleted
	//make a new mrt so it isn't a deep copy (and we can always build off of the original empty one)
	core::chemical::MutableResidueType my_mrt( dummylig_mrt );

	//declare strings to use for the atom name, type, and mm_type
	std::string atom_name = "";
	std::string atom_type_name = "";
	std::string mm_atom_type_name = "";

	//declare vector to hold atom coordinates
	utility::vector1<core::Real> atom_xyz(3,0);

	//declare a vector to hold the names of first 3 atoms so we can apply icoor data (especialy for the first 3 atoms)
	utility::vector1<std::string> first_three_atoms;

	//counter to determine the first 3 atoms so we can assign icoor_data to the first 3 atoms after we reach them (can't do beforehand, as we do not know atoms before we reach them)
	core::Size atom_counter = 0;

	//counter to count the total atom number for unique name tracking
	core::Size total_atom_counter = 0;

	//create pose for matrix
	pose::Pose matrix_pose;

	//iterate through the matrix cells to add atoms
	//run through the boundaries and determine the space fill ratio of the sub-area
	for ( core::Size x = 1; x <= xyz_bound[1]; ++x ) {
		for ( core::Size y = 1; y <= xyz_bound[2]; ++y ) {
			for ( core::Size z = 1; z <= xyz_bound[3]; ++z ) {

				//opt to continue if the atom in question is not wanted for the export per flags
				//use_empty_space_in_space_fill_pdb or use_whole_matrix_in_space_fill_pdb
				//cases where the atom is indicated to be outside of the sub area
				if ( option[ OptionKeys::motifs::use_whole_matrix_in_space_fill_pdb ]() == false && (space_fill_matrix[x][y][z] == 0 || space_fill_matrix[x][y][z] == 1 || space_fill_matrix[x][y][z] == 5 || space_fill_matrix[x][y][z] == 9 ) ) {
					continue;
				}
				//cases where we do not want atoms to represent the empty space
				//values that are even correspond to empty space (currently just 0 and 2)
				if ( option[ OptionKeys::motifs::use_empty_space_in_space_fill_pdb ]() == false && space_fill_matrix[x][y][z] % 2 == 0 ) {
					continue;
				}

				//derive unique atom name
				atom_name = base_10_to_base_62(total_atom_counter);

				//adjust the atom type based on the value in the matrix cells
				/*
				0 = empty and out of sub area, carbon, black
				1 = protein and out of sub area, fluorine, icy blue
				2 = empty and in sub area, oxygen, red
				3 = protein and in sub area, nitrogen, blue
				4 = do not use, keep even for unoccupied space
				5 = ligand and out of sub area, sulphur, yellow
				6 = do not use, keep even for unoccupied space
				7 = ligand and in sub area, chlorine, green
				8 = do not use, keep even for unoccupied space
				9 = ligand and protein out of sub area, phosphorous, orange
				10 = do not use, keep even for unoccupied space
				11 = ligand and protein in sub area, iodine, purple
				*/
				//carbon
				if ( space_fill_matrix[x][y][z] == 0 ) {
					atom_type_name = "aroC";
					mm_atom_type_name = "C";
				}
				//hydrogen
				//hydrogen may be a problem, try fluorine instead (icy blue)
				if ( space_fill_matrix[x][y][z] == 1 ) {
					//atom_type_name = "Haro";
					//mm_atom_type_name = "H";
					atom_type_name = "F";
					mm_atom_type_name = "F1";
				}
				//oxygen
				if ( space_fill_matrix[x][y][z] == 2 ) {
					atom_type_name = "OH";
					mm_atom_type_name = "O";
				}
				//nitrogen
				if ( space_fill_matrix[x][y][z] == 3 ) {
					atom_type_name = "NH2O";
					mm_atom_type_name = "N";
				}
				//sulphur
				if ( space_fill_matrix[x][y][z] == 5 ) {
					atom_type_name = "S";
					mm_atom_type_name = "S";
				}
				//chlorine
				if ( space_fill_matrix[x][y][z] == 7 ) {
					atom_type_name = "Cl";
					mm_atom_type_name = "CL";
				}
				//phosphorous
				if ( space_fill_matrix[x][y][z] == 9 ) {
					atom_type_name = "Pha";
					mm_atom_type_name = "P";
				}
				//iodine
				if ( space_fill_matrix[x][y][z] == 11 ) {
					atom_type_name = "I";
					mm_atom_type_name = "I";
				}
				//if true (occupied), nitrogen
				//make the atom
				//for simplicity sake, charge is 0
				my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);

				//derive the coordinates of the atom
				//divide matrix coordinates by resolution scalar and then subtract corresponding shift
				atom_xyz[1] = (static_cast<core::Real>(x)/resolution_increase_factor) - xyz_shift[1];
				atom_xyz[2] = (static_cast<core::Real>(y)/resolution_increase_factor) - xyz_shift[2];
				atom_xyz[3] = (static_cast<core::Real>(z)/resolution_increase_factor) - xyz_shift[3];

				//adjust the coordinates of the atom in the mrt
				//TEST THIS: looks like I can use MutableResidueType::set_ideal_xyz to set the xyz (and hopefully not have to deal with icoor); use function that takes in the atom name string (it just calls the overload the uses the vertex anyway)
				my_mrt.set_ideal_xyz(atom_name, Vector(atom_xyz[1],atom_xyz[2],atom_xyz[3]) );

				//handle adding the icoor data for the atom
				//skip if the atom name is "" (which will happen if you only look at the sub-matrix)

				if ( atom_name != "" ) {
					++atom_counter;
					++total_atom_counter;

					//if atom counter is <= 3, add the atom name to the first_three_atoms vector
					//atom count must be >3 so we know the first 3 atoms and can retroactively add icoor for those first 3 atoms
					if ( atom_counter <= 3 ) {
						first_three_atoms.push_back(atom_name);

						//if the count is 3, we can retroactively add icoor data for the first 3 atoms now
						if ( atom_counter == 3 ) {

							my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);

							my_mrt.set_icoor(first_three_atoms[2],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
							//flip order of stub2 and stub3 for 3rd atom (looking at some params files, it seems like they do this)
							my_mrt.set_icoor(first_three_atoms[3],0,0,0,first_three_atoms[1],first_three_atoms[3],first_three_atoms[2],false);

							//add bonds to connect these 3
							my_mrt.add_bond(first_three_atoms[1],first_three_atoms[2]);
							my_mrt.add_bond(first_three_atoms[2],first_three_atoms[3]);
						}
					} else {
						//add icoor data for later atoms, using the first 3 atoms as the same stub atoms each time
						my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);

						//bond everything else to the 3rd atom and hope things dont get wonky
						//if they do, I can try another means to attach atoms to each other looking at previous atoms
						my_mrt.add_bond(first_three_atoms[3],atom_name);

					}

					//cut off when my_mrt becomes 100 atoms and start a new residue so that the conversion for an rt doesn't take too long
					//there might need to be some handling at the end if my_mrt only has 1-2 atoms
					//handling would probably just be to add fake atoms in the place of the first atom to fill out to 3 for the icoor data
					if ( atom_counter == 100 ) {
						//convert my_mrt to a data type that can be added to a pose

						//assign internal coordinates
						my_mrt.assign_internal_coordinates();

						//convert to residuetype
						ResidueTypeCOP my_rt(ResidueType::make(my_mrt));

						//convert to residue
						core::conformation::ResidueOP my_res( core::conformation::ResidueFactory::create_residue(*my_rt));

						//append residue to pose
						matrix_pose.append_residue_by_jump( *my_res, 1 );

						//reset atom counter
						atom_counter = 0;

						//pop back the entries to first_three_atoms so we can start again
						first_three_atoms.pop_back();
						first_three_atoms.pop_back();
						first_three_atoms.pop_back();

						//reset my_mrt
						my_mrt = dummylig_mrt;
					}

					//set name back to "" so that we don't keep going after we finish the sub area
					atom_name = "";

				}



			}
		}
	}

	
	ms_tr.Trace << "Made all atoms for this matrix" << std::endl;
	

	//handle rare case where the my_mrt that reaches this point only has 1-2 atoms and would have no icoor data
	//if no atoms, move to dumping the pdb
	if ( my_mrt.natoms() == 0 ) {
		
		ms_tr.Trace << "Making viualization pdb " << matrix_pdb_name << std::endl;
		
		matrix_pose.dump_pdb(matrix_pdb_name);

		return matrix_pose;
	} else if ( my_mrt.natoms() == 1 ) {
		//if only 1 atom, make 2 "copies" (atom in the same position as original) of the atom and assign icoor data
		//make name of atom
		atom_name = base_10_to_base_62(total_atom_counter + 1);
		//add atom to mrt
		//use same types as what was used last
		my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);
		//add atom to the first three atoms list for tracking
		first_three_atoms.push_back(atom_name);

		//make name of atom
		atom_name = base_10_to_base_62(total_atom_counter + 2);
		//add atom to mrt
		//use same types as what was used last
		my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);
		//add atom to the first three atoms list for tracking
		first_three_atoms.push_back(atom_name);

		//set icoor data for these atoms
		my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[2],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[3],0,0,0,first_three_atoms[1],first_three_atoms[3],first_three_atoms[2],false);

		//add bonds to connect these 3
		my_mrt.add_bond(first_three_atoms[1],first_three_atoms[2]);
		my_mrt.add_bond(first_three_atoms[2],first_three_atoms[3]);

	} else if ( my_mrt.natoms() == 2 ) {
		//if 2 atoms, make a 3rd atom that is a "copy" of the first atom
		//make name of atom
		atom_name = base_10_to_base_62(total_atom_counter + 1);
		//add atom to mrt
		//use same types as what was used last
		my_mrt.add_atom(atom_name, atom_type_name, mm_atom_type_name, 0);
		//add atom to the first three atoms list for tracking
		first_three_atoms.push_back(atom_name);

		//set icoor data for these atoms
		my_mrt.set_icoor(first_three_atoms[1],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[2],0,0,0,first_three_atoms[1],first_three_atoms[2],first_three_atoms[3],false);
		my_mrt.set_icoor(first_three_atoms[3],0,0,0,first_three_atoms[1],first_three_atoms[3],first_three_atoms[2],false);

		//add bonds to connect the 3 atoms
		my_mrt.add_bond(first_three_atoms[1],first_three_atoms[2]);
		my_mrt.add_bond(first_three_atoms[2],first_three_atoms[3]);
	}
	//now assign internal coordina

	//try to use assign_internal_coordinates to see if that helps fix icoor
	my_mrt.assign_internal_coordinates();

	//convert my_mrt to a data type that can be added to a pose
	//convert to residuetype
	ResidueTypeCOP my_rt(ResidueType::make(my_mrt));

	//convert to residue
	core::conformation::ResidueOP my_res( core::conformation::ResidueFactory::create_residue(*my_rt));



	//append residue to pose
	matrix_pose.append_residue_by_jump( *my_res, 1 );

	//export pose using custom name
	//custom name based on prefix (i.e. motif used, "empty", etc.), occupied ratios


	
	ms_tr.Trace << "Making viualization pdb " << matrix_pdb_name << std::endl;
	
	matrix_pose.dump_pdb(matrix_pdb_name);

	return matrix_pose;
}


// @brief function to be used to convert a base 10 number to base 62 (as a string with characters represented by base_62_cipher_)
//used in export_space_fill_matrix_as_C_H_O_N_pdb to assign a unique name to an atom (due to limitations in atom icoor data, an atom name can be no longer than 4 characters)
std::string LigandDiscoverySearch::base_10_to_base_62(core::Size starting_num)
{
	//bool to indicate whether to keep processing the number
	//stop processing once we fully build the string, which occurs when dividing the current number by 62 is <1 (integer would be 0)
	bool keep_processing = true;

	//string to hold the base 62 representation of the number
	//characters represented by base_62_cipher_ vector in the .hh file
	std::string base_62_number = "";

	core::Size curr_num = starting_num;

	while ( keep_processing )
	{
		//take the quotient of the current number by 62
		core::Size quotient = curr_num / 62;

		//derive the modulus of the current number by 62
		core::Size mod = curr_num % 62;

		//use the mod value to get the corresponding character from the cipher and append to the string
		//character appends to the front of the string
		base_62_number.insert(base_62_number.begin(),utility::code_to_6bit(mod + 1));

		//if the quotient is under 62, get the number from the cipher, otherwise we have to repeat the processing operation
		if ( quotient < 62 ) {
			// use the quotient to add to the number, unless the quotient is 0 (no need to add a placeholder 0)
			if ( quotient != 0 ) {
				base_62_number.insert(base_62_number.begin(),utility::code_to_6bit(quotient + 1));
			}

			//we have now fully derived the base 62 number and can stop processing
			keep_processing = false;
		} else {
			//set the quotient to the current number and we continue the operation off of it
			curr_num = quotient;
		}
	}

	return base_62_number;
}

// @brief function to hash out an input motif library into a standard map of motifCOPs
//inputs are initial motif library and map that is to be filled out
//map keys are tuples of 7 strings, which is the residue involved in the motif and then the names of the atoms involved (3 atoms on both sides of motif; we don't care about ligand name in key)
//
void LigandDiscoverySearch::hash_motif_library_into_map(protocols::motifs::MotifCOPs & input_library, std::map<motif_atoms, protocols::motifs::MotifCOPs> & mymap)
{
	//iterate over the input library
	for ( auto motifcop : input_library ) {
		//collect residue name from motifcop
		std::string motif_residue_name(motifcop->restype_name1());
		//declare tuple key
		motif_atoms key_tuple(motifcop->restype_name1(),motifcop->res1_atom1_name(),motifcop->res1_atom2_name(),motifcop->res1_atom3_name(),motifcop->res2_atom1_name(),motifcop->res2_atom2_name(),motifcop->res2_atom3_name());

		//push back motif to motifcops at key address
		mymap[key_tuple].push_back(motifcop);

	}
}