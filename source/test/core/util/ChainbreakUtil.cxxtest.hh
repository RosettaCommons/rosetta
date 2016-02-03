// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/util/ChainbreakUtil.cxxtest.hh
/// @author Christopher Miles (cmiles@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Utility headers
#include <numeric/xyzVector.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>

using core::pose::Pose;

class ChainbreakUtilTest : public CxxTest::TestSuite {
 public:
  Pose pose_;

  void setUp() {
    core_init();
    core::import_pose::pose_from_file(pose_, "core/kinematics/test.pdb", core::import_pose::PDB_file);
    core::util::switch_to_residue_type_set(pose_, core::chemical::CENTROID);
  }

  void test_has_chainbreak() {
    using core::id::NamedAtomID;

    NamedAtomID id_n("N", 1);
    NamedAtomID id_c("C", 1);
    NamedAtomID id_ca("CA", 1);

    // Relocate backbone atoms of 1st residue
    Pose copy(pose_);
    copy.set_xyz(id_n, copy.xyz(id_n) * 2);
    copy.set_xyz(id_c, copy.xyz(id_c) * 2);
    copy.set_xyz(id_ca, copy.xyz(id_ca) * 2);
  }
};
