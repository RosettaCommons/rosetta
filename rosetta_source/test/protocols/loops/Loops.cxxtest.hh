// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/Loops.cxxtest.hh
/// @brief test suite for protocols/loops/Loops
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <numeric/xyzVector.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

namespace {

#define TOLERANCE 0.00001

using protocols::loops::Loop;
using protocols::loops::Loops;

class LoopsTest : public CxxTest::TestSuite {
 public:
  core::pose::PoseOP pose_;
  Loops loops_;

  void setUp() {
    core_init();
    pose_ = core::import_pose::pose_from_pdb("protocols/nonlocal/2GB3.pdb");
    loops_.push_back(Loop(1, 5));
    loops_.push_back(Loop(10, 15));
    loops_.push_back(Loop(20, 25));
  }

  void test_center_of_mass() {
    using core::Real;
    using numeric::xyzVector;

    // Expected
    Real x = -1.6075882353;
    Real y = -3.0146470588;
    Real z = -0.0540588235;

    // Actual
    xyzVector<Real> center;
    loops_.center_of_mass(*pose_, &center);

    TS_ASSERT_DELTA(x, center.x(), TOLERANCE);
    TS_ASSERT_DELTA(y, center.y(), TOLERANCE);
    TS_ASSERT_DELTA(z, center.z(), TOLERANCE);
  }
};

}  // anonymous namespace
