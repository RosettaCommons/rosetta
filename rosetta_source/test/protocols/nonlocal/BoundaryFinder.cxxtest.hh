// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BoundaryFinder.cxxtest.hh
/// @brief test suite for protocols/nonlocal/BoundaryFinder
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <protocols/nonlocal/BoundaryFinder.hh>

namespace {
using core::Size;
using core::kinematics::FoldTree;
using core::pose::Pose;

class BoundaryTest : public CxxTest::TestSuite {
 public:
	Pose pose_;

	void setUp() {
    core_init();
    core::import_pose::pose_from_pdb(pose_, "protocols/nonlocal/2GB3.pdb");
	}

	void tearDown() {}

	void test_simple_boundaries() {
		Size position = 4;  // arbitrary residue position
		Size lower, upper;
		protocols::nonlocal::BoundaryFinder::boundaries(pose_.fold_tree(),
                                                                position,
                                                                &lower,
                                                                &upper);

		TS_ASSERT_EQUALS(lower, 1);
		TS_ASSERT_EQUALS(upper, pose_.total_residue());
	}
};
}  // anonymous namespace
