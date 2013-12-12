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
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/rms_util.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/exit.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <algorithm>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {

void setup_fold_tree(Pose & pose, Loop const & loop) {
	using core::kinematics::FoldTree;
	using core::kinematics::Edge;
	using core::conformation::Residue;
	using core::chemical::CUTPOINT_LOWER;
	using core::chemical::CUTPOINT_UPPER;
	using core::pose::add_variant_type_to_pose_residue;

	// Create and install the fold tree.
	
	FoldTree tree;

	Size start = loop.start() - 1;
	Size stop = loop.stop() + 1;
	Size length = pose.total_residue();

	if (start == 1 || stop == length) {
		tree.simple_tree(length);
	} else {
		tree.add_edge(1, start, Edge::PEPTIDE);
		tree.add_edge(start, stop, Edge::PEPTIDE);
		tree.add_edge(stop + 1, length, Edge::PEPTIDE);
		tree.add_edge(start, stop + 1, 1);
	}

	Size root = (start == 1) ? length : 1;
	tree.reorder(root);

	pose.fold_tree(tree);

	// Set cutpoint variants for correct chain-break scoring.

	if ( stop != length ) {
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

