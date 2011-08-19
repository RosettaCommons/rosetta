// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jumping/Dssp.cxxtest.hh
/// @brief  test suite for Dssp secondary structure calculations.
/// @author James Thompson

#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>

#include <protocols/jumping/util.hh>

#include <string>

//Auto Headers
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/FragCache.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedStubID.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <utility/stream_util.hh>
#include <protocols/Protocol.fwd.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>



static basic::Tracer TR("test.protocols.jumping.DsspTests");

using namespace core;

class DsspTests : public CxxTest::TestSuite
{
public:
	DsspTests() {}

	// shared data
	pose::Pose pose_1ten, pose_test_in;

	// Shared initialization goes here.
	void setUp() {
		core_init();

		pose_test_in = create_test_in_pdb_pose(); // alpha
		pose_1ten    = create_1ten_pdb_pose(); // beta
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_dssp() {
		using core::Size;
		using namespace protocols::jumping;
		//Dssp dssp( pose_1ten );
		//dssp.insert_ss_into_pose( pose_1ten );
		assign_ss_dssp( pose_1ten );
		TS_ASSERT_EQUALS(
			pose_1ten.secstruct(),
			"LLLLEEEEEELLLLLLEEEEEELLLLLLLEEEEEEEELLLLLLLEEEEEELLLLEEEELLLLLLLEEEEEEEEEELLEELLLEEEEEEL"
		);

		//Dssp test_in_dssp( pose_test_in );
		//test_in_dssp.insert_ss_into_pose( pose_test_in );
		assign_ss_dssp( pose_test_in );
		TS_ASSERT_EQUALS(
			pose_test_in.secstruct(),
			"LHHHHHHHHHHHHLLLLLLLLLHHHHHHLLLLHHHHHHHHHHHHLLLHHHHHHHHHHHHHHHHHHHLLLLHHHHHHHHLLLLHHHHHHHHHHHHLLLHHHHHLLLLLLLLLLLLLL"
		);
	}

}; // DsspTests
