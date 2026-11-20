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

// Core
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>

// Motifs
#include <protocols/motifs/Motif.fwd.hh>

// Utility Headers
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>

// c++ headers
#include <map>
#include <string>

namespace protocols {
namespace motifs {

class MotifCloud
{
public:

	/// @brief Integer pose-space voxel coordinate; each voxel spans one angstrom per axis.
	using VoxelIndex = utility::fixedsizearray1<int, 3>;
	struct VoxelIndexLess {
		bool operator()( VoxelIndex const & lhs, VoxelIndex const & rhs ) const {
			for ( core::Size axis = 1; axis <= 3; ++axis ) {
				if ( lhs[axis] != rhs[axis] ) {
					return lhs[axis] < rhs[axis];
				}
			}
			return false;
		}
	};
	/// @brief Counts keyed by 1-based fa_standard atom type index.
	using AtomTypeCounts = std::map<core::Size, core::Size>;
	using SparseMotifMatrix = std::map<VoxelIndex, AtomTypeCounts, VoxelIndexLess>;
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
