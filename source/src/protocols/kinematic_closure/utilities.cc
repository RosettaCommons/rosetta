// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Source file for the solution picker functions.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/utilities.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/VariantType.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

// Utility headers
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <algorithm>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {

/// @details It is not absolutely necessary to call this function before 
/// invoking kinematic closure, but doing so has some advantages compared to 
/// using a simple one-edge fold tree or a standard loop modeling fold tree.  
/// For conciseness, I'll refer to the fold trees created by this function as 
/// 'KIC trees' and I'll refer to the fold trees created by the two alternative 
/// methods as 'simple trees' and 'loop trees' respectively.  KIC and loop 
/// trees are both more efficient than simple trees, because fewer coordinates 
/// have to be recalculated after each move.  The former trees also prevent 
/// loop moves from propagating throughout the whole structure, but in practice 
/// this is a non-issue because kinematic closure is far more precise than the 
/// PDB file format, which only holds three decimal places.  The advantage of a 
/// KIC trees over loop trees is that KIC trees don't create a cutpoint in the 
/// middle of the loop.  Since degrees of freedom that cross the cutpoint can't 
/// be read (which is weird; I might be doing something wrong), KIC isn't able 
/// to read default values for those DOFs from the loop.  This doesn't matter 
/// if all degrees of freedom (including bond angles and lengths) are perturbed 
/// explicitly, but otherwise this could lead to confusing bugs.
///
/// This method is kept in the kinematic_closure namespace rather than the 
/// loop_modeling namespace because the details of the KIC algorithm make 
/// unique demands on the topology of the fold tree.  However, this leads to a 
/// couple of problems:
/// 
/// 1. In some cases the fold tree will get built way more often than 
///    necessary.  (And doing so is slow.)  Consider the case where multiple 
///    loops are being sampled.  The loop_modeling framework knows this and can 
///    build the fold tree just once.  The KIC framework doesn't, and will 
///    build a new single-loop fold tree every time the loops are switched.
///
/// 2. The KIC framework doesn't have a way to restore the original fold tree.  
///    This could be lead to counter-intuitive and annoying bugs in code that 
///    uses KIC as part of a larger protocol and doesn't expect the fold tree 
///    to be clobbered every time KIC is called.
///
/// The punchline is that we need a smarter way to handle fold trees in KIC.  
/// This will probably involve an option in KicMover to specify whether or not 
/// the fold tree needs to be set.  This option would be on by default, but 
/// disabled by the loop modeling protocol.

/*
void setup_fold_tree(Pose & pose, Loop const & loop) {
	using core::kinematics::FoldTree;
	using protocols::loops::fold_tree_from_loops;
	using protocols::loops::Loops;

	Loops loops; loops.add_loop(loop);
	FoldTree tree; fold_tree_from_loops(pose, loops, tree, true);
	pose.fold_tree(tree);
}
*/

void setup_fold_tree(Pose & pose, Loop const & loop) {
	using core::kinematics::FoldTree;
	using core::kinematics::Edge;
	using core::conformation::Residue;
	using core::chemical::CUTPOINT_LOWER;
	using core::chemical::CUTPOINT_UPPER;
	using core::pose::add_variant_type_to_pose_residue;
	using core::conformation::symmetry::get_asymm_unit_fold_tree;

	FoldTree tree;
	FoldTree const original_tree = get_asymm_unit_fold_tree(pose.conformation());

	Size num_residues = original_tree.nres();
	Size num_ligands = 0;

	// We need to know if there are any ligands in the pose, because ligands 
	// cannot be linked into the fold tree by peptide edges.  I'm assuming that 
	// rosetta always loads ligands onto the end of the residue list.

	while (not pose.residue(num_residues).is_protein()) {
		num_residues -= 1;
		num_ligands += 1;
	}

	// Configure the loop section of the fold tree so that the entire loop is 
	// represented by one edge.  Not having any cutpoints in the loop makes it 
	// possible to extract default values for all torsion angles.  This isn't 
	// necessary as long as you're using a perturber that explicitly sets all 
	// DOFs (the most common use case), but it's nice to not have any surprises. 

	Size start = loop.start() - 1;
	Size stop = loop.stop() + 1;
	Size jump_id = 1;

	if (start == 1 || stop == num_residues) {
		tree.simple_tree(num_residues);
	} else {
		tree.add_edge(1, start, Edge::PEPTIDE);
		tree.add_edge(start, stop, Edge::PEPTIDE);
		tree.add_edge(stop + 1, num_residues, Edge::PEPTIDE);
		tree.add_edge(start, stop + 1, jump_id++);
	}

	// If there are any ligands in this pose, attach them to the main protein 
	// chain via jumps.

	for (Size i = num_residues + 1; i <= num_residues + num_ligands; i++) {
		tree.add_edge(1, i, jump_id++);
	}

	// Finish building the fold tree and assign it to the given pose.

	tree.reorder(1);
	pose.fold_tree(tree);

	// Set cutpoint variants for correct chain-break scoring.

	if ( stop != num_residues ) {
		Residue const & lower_residue = pose.residue(stop);
		Residue const & upper_residue = pose.residue(stop + 1);

		if (lower_residue.has_variant_type(CUTPOINT_LOWER) == false)
			add_variant_type_to_pose_residue(pose, CUTPOINT_LOWER, stop);
		if (upper_residue.has_variant_type(CUTPOINT_UPPER) == false)
			add_variant_type_to_pose_residue(pose, CUTPOINT_UPPER, stop + 1);
	}
}

} // namespace kinematic_closure
} // namespace protocols

