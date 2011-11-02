// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopsUtil.cxxtest.hh
/// @brief test suite for protocols/loops/util
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

//Auto Headers
#include <utility/vector1.hh>


namespace {

using core::Size;
using core::pose::Pose;
using core::pose::PoseOP;
using protocols::loops::Loop;
using protocols::loops::Loops;

class LoopsUtilTest : public CxxTest::TestSuite {

private:
  bool equal_torsions(const Pose& p1, const Pose& p2) {
    if (p1.total_residue() != p2.total_residue())
      return false;

    for (Size i = 1; i <= p1.total_residue(); ++i) {
      if (p1.phi(i) != p2.phi(i)) {
        return false;
      } else if (p1.psi(i) != p2.psi(i)) {
        return false;
      } else if (p1.omega(i) != p2.omega(i)) {
        return false;
      }
    }

    return true;
  }

public:
  void setUp() {
    protocols_init();
  }

  void testSafeExtendLoopsAndIdealize() {
    PoseOP pose = core::import_pose::pose_from_pdb("protocols/loops/2GB3.pdb");
    Pose other = *pose;

    Loops loops;
    protocols::loops::safe_set_extended_torsions_and_idealize_loops(loops, &other);
    TS_ASSERT(equal_torsions(*pose, other));

    loops.push_back(Loop(1,5));
    protocols::loops::safe_set_extended_torsions_and_idealize_loops(loops, &other);
    TS_ASSERT(!equal_torsions(*pose, other));
  }
};

}  // anonymous namespace
