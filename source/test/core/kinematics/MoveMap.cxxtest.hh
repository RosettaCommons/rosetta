// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/MoveMap.cxxtest.hh
/// @brief  test suite for core::kinematics::MoveMap.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/simple_moves/BackboneMover.hh>

// Utility headers
#include <basic/Tracer.hh>

// C/C++
#include <iostream>
#include <vector>

//Auto Headers
#include <utility/vector1.hh>


using basic::Tracer;
using core::kinematics::MoveMap;
using core::kinematics::MoveMapOP;
using core::pose::Pose;
using std::pair;
using std::vector;

namespace {
Tracer TR("core.kinematics.MoveMap.cxxtest");

class MoveMapTest : public CxxTest::TestSuite {
 public:
  Pose pose_;

  void setUp() {
    core_init();
    core::import_pose::pose_from_pdb(pose_, "core/kinematics/test.pdb");
  }

  void tearDown() {}

  // Ensure that backbone torsions that we wish to remain fixed are not
  // altered through fragment insertion operations
  void test_backbone_range_protection() {
		using core::Real;
		using core::Size;
		using protocols::simple_moves::SmallMover;
		using std::endl;

		vector<Real> before;
		for (Size i = 11; i <= 19; ++i) {
			before.push_back(pose_.phi(i));
			before.push_back(pose_.psi(i));
			before.push_back(pose_.omega(i));
		}

		// identify portions of the backbone whose torsions should remain unchanged
		vector<pair<Size, Size> > offlimits;
		offlimits.push_back(std::make_pair(11, 19));

		MoveMapOP mmap( new MoveMap() );
		mmap->set_ranges_unmodifiable(offlimits);

		// create a simple mover and make some modifications
		Pose modified_pose(pose_);
		protocols::simple_moves::SmallMover mover(mmap, 10, 200);
		mover.apply(modified_pose);

		// make sure that some moves have occurred
		Real rmsd = core::scoring::CA_rmsd(pose_, modified_pose);
		TR << "rmsd => " << rmsd << endl;
		TS_ASSERT_DIFFERS(0, rmsd);

		vector<Real> after;
		for (Size i = 11; i <= 19; ++i) {
			after.push_back(modified_pose.phi(i));
			after.push_back(modified_pose.psi(i));
			after.push_back(modified_pose.omega(i));
		}

		// ensure that the protected torsions have not been modified
		TS_ASSERT_EQUALS(before.size(), after.size());
		for (size_t i = 0; i < before.size(); ++i) {
			TR << before[i] << ", " << after[i] << endl;
			TS_ASSERT_EQUALS(before[i], after[i]);
		}
  }
};
}  // anonymous namespace
