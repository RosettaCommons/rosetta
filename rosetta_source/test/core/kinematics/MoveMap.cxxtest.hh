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
// AUTO-REMOVED #include <core/fragment/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
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
// AUTO-REMOVED #include <utility>
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

		MoveMapOP mmap = new MoveMap();
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

	void test_levels_of_specification_1() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		// DOF Level specification
		movemap.set( PHI, true );
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get_chi( 2 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(movemap.get( PHI ));
		TS_ASSERT(movemap.get( DOF_ID( AtomID(1, 1), PHI)));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D)));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_2() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_bb( true );
		TS_ASSERT(movemap.get( BB ));
		TS_ASSERT(movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get_chi( 2 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI)));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D)));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_3() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_bb( 1, true );
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get_chi( 2 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_4() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		utility::vector1<bool> vec_1;
		vec_1.push_back(true);
		movemap.set_bb( vec_1 );
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get_chi( 2 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_5() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_bb_true_range(1, 1);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get_chi( 2 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_6() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_bb( 1 );
		std::vector<std::pair<core::Size, core::Size> > range;
		range.push_back(std::pair<core::Size, core::Size>(1,1));
		movemap.set_ranges_unmodifiable(range);
		TS_ASSERT(movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(movemap.get_bb( 2 ));
		TS_ASSERT(!movemap.get_chi( 2 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_7() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;


		movemap.set_chi( true );
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(movemap.get_chi( 2 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(movemap.get(CHI));
		TS_ASSERT(movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_8() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		utility::vector1<bool> vec_1;
		vec_1.push_back(true);
		movemap.set_chi(vec_1);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_9() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_chi_true_range(1, 1);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_10() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_jump(true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(movemap.get( JUMP ));
		TS_ASSERT(movemap.get_jump( 10 ));
		TS_ASSERT(movemap.get_jump( 10, 11));
		TS_ASSERT(movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();


	}

	void test_levels_of_specification_11() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_jump(10, 11, true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get( JUMP ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(movemap.get_jump( 10, 11));
		TS_ASSERT(movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_12() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_jump(JumpID(10, 11), true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get( JUMP ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(movemap.get_jump( 10, 11));
		TS_ASSERT(movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_13() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;


		movemap.set(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI), true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get( JUMP ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_14() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_jump(true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(movemap.get( JUMP ));
		TS_ASSERT(movemap.get_jump( 10 ));
		TS_ASSERT(movemap.get_jump( 10, 11));
		TS_ASSERT(movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_15() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_jump(10, true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get( JUMP ));
		TS_ASSERT(movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_15_1() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_jump(10, 11, true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get( JUMP ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(movemap.get_jump( 10, 11));
		TS_ASSERT(movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_15_2() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set_jump(JumpID(10, 11), true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get( JUMP ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(movemap.get_jump( 10, 11));
		TS_ASSERT(movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

	void test_levels_of_specification_16() {
		using core::kinematics::MoveMap;
		using namespace core::id;

		MoveMap movemap;

		movemap.set( DOF_ID( AtomID(1, 1), PHI), true);
		TS_ASSERT(!movemap.get( BB ));
		TS_ASSERT(!movemap.get_bb( 1 ));
		TS_ASSERT(!movemap.get( CHI ));
		TS_ASSERT(!movemap.get_chi( 1 ));
		TS_ASSERT(!movemap.get( JUMP ));
		TS_ASSERT(!movemap.get_jump( 10 ));
		TS_ASSERT(!movemap.get_jump( 10, 11));
		TS_ASSERT(!movemap.get_jump(JumpID(10, 11)));
		TS_ASSERT(!movemap.get(CHI));
		TS_ASSERT(!movemap.get(core::kinematics::MoveMap::MoveMapTorsionID(10, CHI)));
		TS_ASSERT(!movemap.get( PHI ));
		TS_ASSERT(movemap.get( DOF_ID( AtomID(1, 1), PHI )));
		TS_ASSERT(!movemap.get( DOF_ID( AtomID(3, 4), D )));
		TS_ASSERT(movemap != *(new MoveMap()));
		TS_ASSERT(movemap == *(movemap.clone()));
		movemap.clear();

	}

};

}  // anonymous namespace
