// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Test the algorithms that drive the KinematicMover.
/// @author Kale Kundert

#ifndef INCLUDED_protocols_kinematic_closure_KinematicClosure_CXXTEST_HH
#define INCLUDED_protocols_kinematic_closure_KinematicClosure_CXXTEST_HH

#define private public
#define protected public

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/protocols/kinematic_closure/TestHelpers.hh>
#include <test/protocols/kinematic_closure/ClosureTests.hh>

// Unit headers
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/utilities.hh>

// Core headers
#include <core/pose/Pose.hh>

// Utility headers
#include <boost/foreach.hpp>
#include <numeric/xyzVector.hh>
#define foreach BOOST_FOREACH

// C++ headers
#include <vector>

using namespace std;
using namespace core;
using protocols::kinematic_closure::ClosureProblem;
using protocols::kinematic_closure::ClosureProblemOP;
using protocols::kinematic_closure::ClosureSolutionCOP;
using protocols::kinematic_closure::SolutionList;

// These test have a couple of shortcomings.  Since I'm not going to fix these 
// things right now, I will at least mention them so they won't cause any 
// surprises later:
//
// 1. There are no tests for cases that can't be closed.  This is because the 
//    closure tests are defined by the atomic coordinates of closed loops.  But 
//    a good way to fix this would be to allow optional internal coordinates to 
//    be specified in the closure problem.  If not given, these coordinates 
//    would simply be taken from the given atomic coordinates.  The mover would 
//    then use these internal coordinates instead of the atomic ones used now.
//
// 2. There are no fold tree tests.

class KinematicClosureTests : public CxxTest::TestSuite {

public:

	void setUp() { core_init(); }

	void test_kinematic_closure() {
		vector<TestHelperOP> test_helpers;

		test_helpers.push_back(new FiveResidueTest);
		test_helpers.push_back(new SixResidueTest);
		test_helpers.push_back(new SevenResidueTest);
		test_helpers.push_back(new FoldTreeTest);
		test_helpers.push_back(new NumericalStabilityTest);

		foreach(TestHelperOP helper, test_helpers) {
			TestMoverOP mover = new TestMover(helper);
			ClosureProblemOP problem = new ClosureProblem();
			SolutionList solutions;

			helper->setup();

			//write_pdb(helper->pose, "closed.pdb");

			problem->perturb(helper->pose, helper->loop, mover);
			problem->solve(solutions);

			helper->test_mover(helper->pose);
			helper->test_closure(solutions);
			helper->test_jacobian(solutions);
			helper->test_pickers(problem, solutions);

			//write_pdb(helper->pose, "open.pdb");
			//write_pdbs(helper->idealized_pose, solutions, "solution");

			problem->restore(helper->pose);
			helper->test_restore(helper->pose);

		}
	}
};

#endif


