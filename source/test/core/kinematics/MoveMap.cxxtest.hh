// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
		core::import_pose::pose_from_file(pose_, "core/kinematics/test.pdb", core::import_pose::PDB_file);
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
		for ( Size i = 11; i <= 19; ++i ) {
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
		for ( Size i = 11; i <= 19; ++i ) {
			after.push_back(modified_pose.phi(i));
			after.push_back(modified_pose.psi(i));
			after.push_back(modified_pose.omega(i));
		}

		// ensure that the protected torsions have not been modified
		TS_ASSERT_EQUALS(before.size(), after.size());
		for ( size_t i = 0; i < before.size(); ++i ) {
			TR << before[i] << ", " << after[i] << endl;
			TS_ASSERT_EQUALS(before[i], after[i]);
		}
	}

	void test_show() {
		MoveMap mm;

		// Do an example of all the settings you can have in the MoveMap
		mm.set_bb(false);
		mm.set_bb(1, true);
		mm.set_chi(true);
		mm.set_chi(2, false);
		mm.set_jump(3, true);
		mm.set_jump(4, false);

		mm.set( core::id::NU, true );
		mm.set( core::id::BRANCH, false );

		mm.set( core::kinematics::MoveMap::MoveMapTorsionID( 5, core::id::NU ), false );
		mm.set( core::kinematics::MoveMap::MoveMapTorsionID( 6, core::id::BRANCH ), true );

		mm.set( core::id::TorsionID(7, core::id::BB, 2), true );
		mm.set( core::id::TorsionID(8, core::id::CHI, 3), false );

		mm.set( core::id::THETA, true );
		mm.set( core::id::DOF_ID( core::id::AtomID(9,10), core::id::D ), true );

		mm.set_jump( core::id::JumpID(11,12) , true );

		std::stringstream out;
		mm.show(out);

		// Feel free to change this string for formatting issues
		// (The comparision is mainly here to make sure all the data gets represented in the output.)
		std::string expected("\n"
			"-------------------------------\n"
			"  resnum     Type  TRUE/FALSE \n"
			"-------------------------------\n"
			" DEFAULT      BB     FALSE\n"
			" DEFAULT      SC      TRUE\n"
			" DEFAULT      NU      TRUE\n"
			" DEFAULT  BRANCH     FALSE\n"
			"     001      BB      TRUE\n"
			"     002      SC     FALSE\n"
			"     005      NU     FALSE\n"
			"     006  BRANCH      TRUE\n"
			"-------------------------------\n"
			" jumpnum     Type  TRUE/FALSE \n"
			"-------------------------------\n"
			" DEFAULT     JUMP    FALSE\n"
			"     003     JUMP     TRUE\n"
			"     004     JUMP    FALSE\n"
			"-------------------------------\n"
			" jumpstart    jumpend  TRUE/FALSE \n"
			"-------------------------------\n"
			"       011        012     TRUE\n"
			"-------------------------------\n"
			"  resnum  atomnum     Type  TRUE/FALSE \n"
			"-------------------------------\n"
			" DEFAULT               PHI    FALSE\n"
			" DEFAULT             THETA     TRUE\n"
			" DEFAULT                 D    FALSE\n"
			" DEFAULT               RB1    FALSE\n"
			" DEFAULT               RB2    FALSE\n"
			" DEFAULT               RB3    FALSE\n"
			" DEFAULT               RB4    FALSE\n"
			" DEFAULT               RB5    FALSE\n"
			" DEFAULT               RB6    FALSE\n"
			"     010      009        D     TRUE\n"
			"-------------------------------\n"
			"  resnum   torsion#     Type  TRUE/FALSE \n"
			"-------------------------------\n"
			"     007        002       BB     TRUE\n"
			"     008        003       SC    FALSE\n"
			"-------------------------------\n"
			"\n");

		TR << "~~~~~~~~FOUND~~~~~~~~~~" << std::endl;
		TR << out.str() << std::endl;
		TR << "~~~~~~~~EXPECTED~~~~~~~" << std::endl;
		TR << expected << std::endl;
		TR << "~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

		TS_ASSERT( out.str() == expected );
	}

};
}  // anonymous namespace
