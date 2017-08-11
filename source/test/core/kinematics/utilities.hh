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

// Core headers
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>

// C++ headers
#include <iostream>
#include <sstream>
#include <iomanip>

void TS_ASSERT_FOLD_TREE_HAS_EDGE( // {{{1
	core::kinematics::FoldTree const & tree, core::Size start, core::Size stop, bool jump=false) {

	for ( core::kinematics::Edge const & edge : tree ) {
		bool start_matches = (edge.start() == start);
		bool stop_matches = (edge.stop() == stop);
		bool jump_matches = (edge.is_jump() == jump);

		// If all these conditions are met, the edge was found and the test passes.
		if ( start_matches && stop_matches && jump_matches ) return;
	}

	// If no matching edge was found, the test fails.
	if ( jump ) std::cout << "No jump '";
	else      std::cout << "No edge '";
	std::cout << std::setw(4) << std::setfill('0') << start;
	std::cout << "--";
	std::cout << std::setw(4) << std::setfill('0') << stop;
	std::cout << "' in the given fold tree:" << std::endl << std::endl;
	tree.show(std::cout);
	TS_FAIL("No match found");
}

void TS_ASSERT_FOLD_TREE_HAS_JUMP( // {{{1
	core::kinematics::FoldTree const & tree, core::Size start, core::Size stop) {

	TS_ASSERT_FOLD_TREE_HAS_EDGE(tree, start, stop, true);
}
// }}}1

#endif


