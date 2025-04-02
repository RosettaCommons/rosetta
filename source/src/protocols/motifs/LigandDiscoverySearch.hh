// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief aginsparg, ipatel, sthyme; this script uses a target receptor, list of motifs, and input ligand(s) to attempt to fit ligands against the receptor against a specified residue index


#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeFinder.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.fwd.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>
#include <protocols/motifs/MotifLibrary.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <protocols/motifs/Motif.fwd.hh>
#include <protocols/motifs/MotifHit.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.fwd.hh>
#include <protocols/motifs/BuildPosition.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/minimization_packing/MinMover.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/PDBInfo.fwd.hh>

#include <core/chemical/GlobalResidueTypeSet.fwd.hh>

#include <protocols/dna/PDBOutput.fwd.hh>

#include <basic/prof.fwd.hh>
#include <basic/Tracer.fwd.hh>



// Utility Headers
#include <utility/io/ozstream.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/excn/Exceptions.fwd.hh>

// c++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <functional>
#include <map>

// option key includes

#include <numeric/xyzVector.fwd.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.fwd.hh>

#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <protocols/ligand_docking/LigandArea.fwd.hh>
#include <protocols/ligand_docking/InterfaceBuilder.fwd.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/ligand_docking/HighResDocker.fwd.hh>
#include <protocols/ligand_docking/LigandDockProtocol.fwd.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/chemical/PoseResidueTypeSet.fwd.hh>

#include <ObjexxFCL/FArray1D.fwd.hh>

#include <protocols/protein_grid/ProteinGrid.fwd.hh>

// Time profiling header
#include <time.h>

//for data type debugging
#include <typeinfo>


////////////////////////////////////////////////////////////////////////////////

class LigandDiscoverySearch
{

public:

	////////////////////////////////////////////////////////////////////////////
	//typedefs to help clean up data types that are really lengthy/convoluted, but not making them classes as reuse purpose may be limited and these are probably pretty specific usages

	// @brief 3D matrix of core::Size values to represent trios of a central atom and 2 adjacent atoms
	// this is used to identify atom trios to iterate over for motifs to potentially align against on the ligand side of the motif
	// this is functionally the same as the SpaceFillMatrix typedef, however this secondary typedef seems useful in avoiding using an object type whose name does not make sense in the given context
	typedef utility::vector1<utility::vector1<utility::vector1<core::Size>>> atom_trios;

	////////////////////////////////////////////////////////////////////////////
	//functions/constructors to set up protocol

	// @brief destructor
	~LigandDiscoverySearch();

	// @brief parameterized constructor to load in motif library, pdb, and ligand library
	LigandDiscoverySearch(core::pose::PoseOP pose_from_PDB, protocols::motifs::MotifCOPs motif_library, utility::vector1<core::conformation::ResidueOP> all_residues, utility::vector1<core::Size> working_position);

	// @brief function to load in a library for the search protocol
	void set_motif_library(protocols::motifs::MotifCOPs motif_library);

	// @brief function to load in a library for the search protocol; read in MotifLibrary object and convert to MotifCOPs
	void set_motif_library(protocols::motifs::MotifLibrary motif_library);

	// @brief return contents of motif_library_
	protocols::motifs::MotifCOPs get_motif_library();

	// @brief function to define the vector of residue indices that we will use for applying motifs for ligand placements
	void set_working_positions(utility::vector1<core::Size> working_position);

	// @brief overload function to define the vector of residue indices that we will use for applying motifs for ligand placements; takes single core::Size
	void set_working_positions(core::Size working_position);

	// @brief function to append an additional vector of working position indices to the existing working_positions_ vector
	void add_working_positions(utility::vector1<core::Size> working_positions);

	// @brief function to append an additional single working position index to the existing working_positions_ vector
	void add_working_positions(core::Size working_positions);

	// @brief return contents of working_positions_
	utility::vector1<core::Size>  get_working_positions();

	// @brief function to load in a pose for the receptor
	void set_working_pose(core::pose::PoseOP pose_from_PDB);

	// @brief function to load in a pose for the receptor, use a regular pose and convert to pointer
	void set_working_pose(core::pose::Pose pose_from_PDB);

	// @brief function to get working_pose_
	core::pose::PoseOP get_working_pose();

	// @brief function to set ligand residues in all_residues_
	void set_all_residues(utility::vector1<core::conformation::ResidueOP> all_residues);

	// @brief function to get all ligand residues
	utility::vector1<core::conformation::ResidueOP> get_all_residues();

	////////////////////////////////////////////////////////////////////////////
	//operation functions

	// @brief main function to run ligand discovery operations
	//needs to have values set for working_pose_, motif_library_, and all_residues_
	//parameter is a string to be a prefix name to use for outputted file names
	core::Size discover(std::string output_prefix);

	// @brief function to get a sub-library of motifs from the main library, based on the residue being used (only get for select residue)
	//function may have use outside discover, so public can use it
	protocols::motifs::MotifCOPs get_motif_sublibrary_by_aa(std::string residue_name);

	// @brief function to push all adjacent atom indices if inputted ligand residue (pointer) into a 3D core::Size vector (atom_trios typedef)
	//this in theory could be useful beyond the scope the discover function, so I will leave this as public
	atom_trios derive_adjacent_atoms_of_ligand(const core::conformation::ResidueOP ligresOP, const core::chemical::AtomTypeSetCOP atset);

	// @brief this function takes in a selected scorefunctionOP and gets the ddg of the selected poseOP (ideally with a placed ligand), and returns the ddg
	// making this a public function, as the usage of this is seems broad enough that it could be used outside discover()
	core::Real get_pose_ddg(core::scoring::ScoreFunctionOP score_fxn, core::pose::PoseOP & my_pose);

private:

	// @brief default constructor
	//will need to use class functions to seed values for input pose, motif library, and ligand library
	//should only use parameterized constructor
	LigandDiscoverySearch();

	// @brief this function is to be called by the constructor(s) to seed initial values to cutoffs that are used for scoring/evaluating metrics of placed ligands in discover() and the functions it calls
	void seed_cutoff_values();

	// @brief prepare score functions for usage in discovery function. Called within discover() and not in a constructor. This probably shouldn't be messed with, so it is kept private
	// This previously just was code in discover(), but it is better to compartmentalize for readability
	void setup_score_functions();

	// @brief function used to make a minipose (focused pose around placed ligand to get quicker scoring of metrics like fa_atr and fa_rep)
	//if returns true, a minipose was successfully made; if returns false, minipose is still empty because no other residues were recruited to it
	bool make_minipose(core::pose::PoseOP & minipose, const core::conformation::ResidueOP ligresOP);

	// @brief function to score the minipose iteratively on fa_rep, fa_atr, and fa_atr+fa_rep
	//this is performed iteratively in order to more quickly filter on select criteria with more and quicker filtering happening at earlier steps
	//returns a boolean based on whether the minipose score was good enough at all steps or not (false means it failed, and the placement will be killed)
	bool score_minipose(const core::pose::PoseOP & minipose, core::Real & fa_rep, core::Real & fa_atr, core::Real & fa_atr_rep_score_before);

	// @brief create a constraint set on the 3 ligand motif atoms (the last residue in working_pose_) to reduce their movement before applying a highresdock
	void add_constraints_to_working_pose(const core::Size trip_atom_1, const core::Size trip_atom_2, const core::Size trip_atom_3, const core::Size working_position, const core::conformation::ResidueOP ligresOP);

	// @brief This function adds a ligand mutableresiduetypeOP to the residue type set for working_pose and original_pose. This occurs on the first instance of this ligand passing enough filters to be used in highresdock
	// the ligand is added to both so that the ligand can be a part of the original_pose for future iterations of placement attempts for the ligand, as well as the working_pose for immediate use
	void add_ligand_to_pose_residuetypeset(const core::chemical::MutableResidueTypeOP lig_mrt);

	// @brief This function creates a HighResDockOP object for using highresdock in discover() to try to optimize the ligand placement. The function takes in a score function OP to use in the HRD
	protocols::ligand_docking::HighResDockerOP make_HighResDockOP_for_discovery(const core::scoring::ScoreFunctionOP my_fxn);

	// @brief This function uses a HighResDockOP object for using highresdock in discover() to try to optimize the ligand placement
	void run_HighResDock_on_working_pose(const protocols::ligand_docking::HighResDockerOP my_HighResDocker);

	// @brief A function to reset the working_pose_ back to its original state, as found in original_pose
	void reset_working_pose();

	//class variables

	// @brief receptor pose to work with
	core::pose::PoseOP working_pose_;
	// @brief original copy of working_pose_, which is used as a copy that working_pose_ can be reverted to after being modified through operations like the highresdock
	// input ligands will be added to the original_pose_ residue type set
	core::pose::Pose original_pose_;
	// @brief motif library (all motifs for all residues)
	protocols::motifs::MotifCOPs motif_library_;
	// @brief motifs library for select residue
	protocols::motifs::MotifCOPs motif_library_for_select_residue_;
	// @brief ligand library
	utility::vector1<core::conformation::ResidueOP> all_residues_;
	// @brief vector to hold list of all indices to investigate/use as anchor residues, used to set value of working_position_
	utility::vector1<core::Size> working_positions_;

	// @brief variable to be used as a cutoff to define the maximum allowed fa_rep score for a placement to be considered
	core::Real fa_rep_cutoff_;

	// @brief variable to be used as a cutoff to define the maximum allowed fa_atr score for a placement to be considered
	core::Real fa_atr_cutoff_;

	// @brief variable to be used as a cutoff to define the maximum allowed combined fa_atr and fa_rep score for a placement to be considered
	core::Real fa_atr_rep_cutoff_;

	// @brief variable to be used as a cutoff to define the maximum allowed score from the whole score function for a placement to be considered
	core::Real whole_fxn_cutoff_;

	// @brief variable to be used as a cutoff to define the maximum allowed ddg score for a placement to be considered
	core::Real ddg_cutoff_;

	// @brief variable to be used as a cutoff to define the minimum number of motif-like contacts that a placed ligand has to be considered
	core::Size min_motifs_cutoff_;

	// @brief variable to be used as a cutoff to define the minimum number of motif-like contacts that a placed ligand has against residues noted as significant to be considered
	core::Size min_sig_motifs_cutoff_;

	// @brief variable to be used as a cutoff to define the minimum number of motif-like contacts that a placed ligand has that are considered "real"
	//real is defined as the motif being within the distance and angle threshold of a motif with the same residue/atoms from the input motif library
	core::Real real_motif_ratio_cutoff_;

	// @brief whole ligand.wts score function to be used with scoring placed ligand poses that pass upstream filters
	core::scoring::ScoreFunctionOP whole_score_fxn_;

	// @brief modified ligand.wts function that only contains fa_atr as a weighted term. Used in quicker preliminary filtering before running whole score function
	core::scoring::ScoreFunctionOP fa_atr_fxn_;

	// @brief modified ligand.wts function that only contains fa_rep as a weighted term. Used in quicker preliminary filtering before running whole score function
	core::scoring::ScoreFunctionOP fa_rep_fxn_;

	// @brief modified ligand.wts function that only contains fa_atr and fa_rep as a weighted term. Used in quicker preliminary filtering before running whole score function
	core::scoring::ScoreFunctionOP fa_atr_rep_fxn_;

	// @brief a copy of either whole_score_fxn_ or fa_atr_rep_fxn_ that is defined by the OptionKeys::motifs::highresdock_with_whole_score_fxn flag and is used in highresdock and scoring operations
	core::scoring::ScoreFunctionOP working_fxn_;

	// @brief a ProteinGrid of the empty working_pose_ to be used for protein-ligand clashing. This does not use the space_fill functionality
	protocols::protein_grid::ProteinGridOP clash_pose_grid_;

	// @brief a ProteinGrid of the empty working_pose_ to be used for space_fill functionality
	//for the space fill analysis, another proteingrid will be cloned of this template in the discover function
	//The clone is because it should be faster to clone to wipe placed ligand data, as opposed to call the wrap and space fill functions again as a form of wiping the object
	protocols::protein_grid::ProteinGridOP sf_pose_grid_;
};
