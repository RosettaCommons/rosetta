// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Sparse voxelized projections of ligand-motif partner anchors around selected pose residues.


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
#include <core/chemical/AtomTypeSet.hh>

#include <core/pose/PDBInfo.fwd.hh>

#include <core/chemical/GlobalResidueTypeSet.fwd.hh>

#include <protocols/dna/PDBOutput.fwd.hh>

#include <basic/prof.fwd.hh>
#include <basic/Tracer.fwd.hh>



// Utility Headers
#include <utility/io/ozstream.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.fwd.hh>

// c++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <functional>
#include <map>
#include <array>

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

namespace protocols {
namespace motifs {

class MotifCloud
{
public:

	/// @brief Integer pose-space voxel coordinate; each voxel spans one angstrom per axis.
	using VoxelIndex = std::array<int, 3>;
	/// @brief Counts keyed by 1-based fa_standard atom type index.
	using AtomTypeCounts = std::map<core::Size, core::Size>;
	using SparseMotifMatrix = std::map<VoxelIndex, AtomTypeCounts>;
	using PerResidueMotifMatrices = std::map<core::Size, SparseMotifMatrix>;

	// destructor
	~MotifCloud();

	//constructors:

	// @brief create motif projections from selected pose residue indices; an empty selection means all residues
	// Projected coordinates and the aggregate use the shared pose XYZ frame.
	MotifCloud(
		core::pose::PoseOP in_pose,
		protocols::motifs::MotifCOPs motif_library,
		utility::vector1<core::Size> const & source_residue_indices = utility::vector1<core::Size>()
	);

	SparseMotifMatrix const & aggregate_matrix() const { return aggregate_motif_matrix_; }
	PerResidueMotifMatrices const & per_residue_matrices() const { return per_residue_matrices_; }
	core::chemical::AtomTypeSetCOP atom_type_set() const { return atom_type_set_; }
	void write_visualization_pdb(SparseMotifMatrix const & matrix, std::string const & filename) const;

private:

	//Functions==========================================================================================

	// @brief seed the atom_type_set_ with the fa_standard library
	void generate_fa_standard_atom_type_set();

	// @brief project motifs from residues in the working pose
	void project_motifs_from_residues(utility::vector1<core::Size> const & source_residue_indices);

	//Member Variables===================================================================================

	/// @brief sparse aggregate projected anchor counts keyed by pose-space voxel and fa_standard atom-type index
	SparseMotifMatrix aggregate_motif_matrix_;

	/// @brief sparse projected anchor counts for each source residue, keyed by its 1-based pose index
	PerResidueMotifMatrices per_residue_matrices_;

	// @brief the pose whose shared XYZ frame is used for voxel coordinates
	core::pose::PoseOP working_pose_;

	// @brief the full input motifs list provided by the user
	protocols::motifs::MotifCOPs motif_library_;

	// @brief a map of smaller motif libraries separated by protein residue for faster processing at the cost of memory overhead
	//will be filled out when iterating over residues to project motifs from, and populated using the motif_utils get_motif_sublibrary_by_aa() function. Keys will be strings for residue 3 letter codes
	std::map< std::string, protocols::motifs::MotifCOPs > motif_library_map_by_residue_;

	// @brief fa_standard atom type set used to interpret atom type indices in voxel counts
	core::chemical::AtomTypeSetCOP atom_type_set_;

};

} // namespace motifs
} // namespace protocols