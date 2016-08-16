// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/util.hh
/// @brief  Kinematics utility functions
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_util_hh
#define INCLUDED_core_kinematics_util_hh


// Package headers
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/types.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>

#include <core/types.hh>

#include <string>
#include <map>


namespace core {
namespace kinematics {

typedef utility::vector1< utility::vector1< Size > > Links;

/// @brief creat an atom and add it to the residue atom-tree based on information stored in links.
tree::AtomOP
add_atom(
	int const atomno,
	int const seqpos,
	utility::vector1< utility::vector1< Size > > const & links,
	AtomPointer1D & atom_ptr,
	bool const add_jump_atom
);

int
pick_loopy_cutpoint(
	Size const n_res,
	ObjexxFCL::FArray1D_float const & cut_bias_sum
);

tree::AtomOP
setup_backrub_atom_tree(
	utility::vector1< id::AtomID > mainchain, // make our own local copy
	id::AtomID const & downstream_id, // mainchain child of last mainchain atom
	AtomPointer2D const & old_atom_pointer,
	utility::vector1< std::pair< Size, Size > > const & edges,
	Size const first_new_pseudo_residue
);

/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
void
simple_visualize_fold_tree( FoldTree const & fold_tree, std::ostream& out );

/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
///                                   **********xxxxxxxxxxxxx************************************
void
simple_visualize_fold_tree_and_movemap( FoldTree const & fold_tree,  MoveMap const& mm, std::ostream& out );

/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
///                                   **********xxxxxxxxxxxxx************************************
void
simple_visualize_fold_tree_and_movemap_bb_chi( FoldTree const & fold_tree, MoveMap const& mm, std::ostream& out );

/// @brief linearizes (or defoliates, if you prefer) a FoldTree.  "default" FoldTrees produced by the PDB reader have all chains (peptide edges) starting from jumps relative to residue 1.  This code modifies the tree to instead have all the jumps be relative to the preceding edge.  It is not tested with ligands and will not work with "functional" jumps.  From A to B:
///A:FOLD_TREE  EDGE 1 78 -1  EDGE 1 79 1   EDGE 79 454 -1  EDGE 1 455 2    EDGE 455 540 -1  EDGE 1 541 3    EDGE 541 697 -1
///B:FOLD_TREE  EDGE 1 78 -1  EDGE 78 79 1  EDGE 79 454 -1  EDGE 454 455 2  EDGE 455 540 -1  EDGE 540 541 3  EDGE 541 697 -1
core::kinematics::FoldTree
linearize_fold_tree( core::kinematics::FoldTree const & tree );

/// sheffler
std::string
visualize_fold_tree(
	FoldTree const & fold_tree
);

std::string
visualize_fold_tree(
	FoldTree const & fold_tree,
	std::map<Size,std::string> const & node_labels_partial
);

std::string
visualize_fold_tree(
	FoldTree const & fold_tree,
	std::map<Size,char> const & mark_jump_to_res
);

std::string
visualize_fold_tree(
	FoldTree const & fold_tree,
	std::map<Size,std::string> const & node_labels_partial,
	std::map<Size,char> const & mark_jump_to_res,
	std::map<Size,Size> const & jump_follows
);

/// @brief remodel a fold tree to account for a large insertion by adding the size of the insert to upstream positions
/// @author Steven Lewis smlewi@gmail.com as a favor for Jared
core::kinematics::FoldTree
remodel_fold_tree_to_account_for_insertion(
	core::kinematics::FoldTree const & input_tree, //return a remodeled version of this tree
	core::Size insert_after, //add insert_size to points after this in primary sequence in the tree
	core::Size insert_size
);

/// @brief Get a vector of residues matching the id from a movemap.
utility::vector1< core::Size >
get_residues_from_movemap_with_id( id::TorsionType query_torsion, MoveMap const & movemap);

/// @brief Get a vector of residues that have any of their BB torsions on - either by way of the full bb or torsion ID setting in Movemap.
utility::vector1< core::Size >
get_residues_from_movemap_bb_any_torsion(MoveMap const & movemap, Size total_resnum);

} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_util_HH
