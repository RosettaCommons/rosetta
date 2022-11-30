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

#include <core/pose/PDBInfo.fwd.hh>

#include <core/chemical/GlobalResidueTypeSet.fwd.hh>

#include <protocols/dna/PDBOutput.fwd.hh>
#include <protocols/dna/util.hh>
#include <protocols/dna/util.hh>
#include <protocols/dna/util.hh>


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


// Time profiling header
#include <time.h>

//for data type debugging
#include <typeinfo>


////////////////////////////////////////////////////////////////////////////////

class LigandDiscoverySearch
{

public:

	////////////////////////////////////////////////////////////////////////////
	//functions/constructors to set up protocol

	// @brief destructor
	~LigandDiscoverySearch();

	// @brief parameterized constructor to load in motif library, pdb, and ligand library
	LigandDiscoverySearch(core::pose::PoseOP pose_from_PDB, protocols::motifs::MotifCOPs motif_library, utility::vector1<core::conformation::ResidueOP> all_residues, core::Size working_position);

	// @brief function to load in a library for the search protocol
	void set_motif_library(protocols::motifs::MotifCOPs motif_library);

	// @brief function to load in a library for the search protocol; read in MotifLibrary object and convert to MotifCOPs
	void set_motif_library(protocols::motifs::MotifLibrary motif_library);

	// @brief return contents of motif_library_
	protocols::motifs::MotifCOPs get_motif_library();

	// @brief function to define the residue index that we will use for applying motifs for ligand placements
	void set_working_position(core::Size working_position);

	// @brief return contents of working_position_
	core::Size get_working_position();

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

private:

	// @brief default constructor
	//will need to use class functions to seed values for input pose, motif library, and ligand library
	//should only use parameterized constructor
	LigandDiscoverySearch();

	// @brief functions to be called in discover() function
	//function to set values for the score functions
	//void seed_score_functions();

	// @brief create protein_representation_matrix_
	//uses working_pose to make the matrix
	void create_protein_representation_matrix(core::Size & x_shift, core::Size & y_shift, core::Size & z_shift, int & x_bound_int, int & y_bound_int, int & z_bound_int);

	// @brief function to run a clash check of the placed ligand against the
	bool ligand_clash_check(core::conformation::ResidueOP ligresOP, core::Size x_shift, core::Size y_shift, core::Size z_shift, int x_bound_int, int y_bound_int, int z_bound_int);

	//class variables

	// @brief receptor pose to work with
	core::pose::PoseOP working_pose_;
	// @brief motif library (all motifs for all residues)
	protocols::motifs::MotifCOPs motif_library_;
	// @brief motifs library for select residue
	protocols::motifs::MotifCOPs motif_library_for_select_residue_;
	// @brief ligand library
	utility::vector1<core::conformation::ResidueOP> all_residues_;
	// @brief residue index of protein pdb to attempt to place ligands off of
	core::Size working_position_;
	// @brief 3D matrix to represent voxelized copy of atoms in pose, used in clash check of placement for quick ruling out of bad placements
	utility::vector1<utility::vector1<utility::vector1<bool>>> protein_representation_matrix_;
};
