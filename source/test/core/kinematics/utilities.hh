// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Header and source file for the TestMover class.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_test_protocols_kinematic_closure_TestHelpers_HH
#define INCLUDED_test_protocols_kinematic_closure_TestHelpers_HH

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/utilities.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>

// Utility headers
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace core;
using core::kinematics::FoldTree;

void TS_ASSERT_FOLD_TREE_HAS_EDGE( // {{{1
		FoldTree const & tree, int start, int stop, bool jump=false) {

	FoldTree::const_iterator it, it_end;

	for (it = tree.begin(), it_end = tree.end(); it != it_end; ++it) {
		bool start_matches = (it->start() == start);
		bool stop_matches = (it->stop() == stop);
		bool jump_matches = (it->is_jump() == jump);

		// If all these conditions are met, the edge was found and the test passes.
		if (start_matches and stop_matches and jump_matches) return;
	}

	// If no matching edge was found, the test fails.
	if (jump) cout << "No jump '";
	else      cout << "No edge '";
	cout << setw(4) << setfill('0') << start;
	cout << "--";
	cout << setw(4) << setfill('0') << stop;
	cout << "' in the given fold tree:" << endl << endl;
	tree.show(cout);
	TS_FAIL("No match found");
}

void TS_ASSERT_FOLD_TREE_HAS_JUMP( // {{{1
		FoldTree const & tree, int start, int stop) {

	TS_ASSERT_FOLD_TREE_HAS_EDGE(tree, start, stop, true);
}
// }}}1

#endif


