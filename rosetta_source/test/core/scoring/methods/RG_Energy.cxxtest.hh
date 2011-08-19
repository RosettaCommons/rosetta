// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   RG_Energy.cxxtest.hh
/// @brief  test suite for RG_Energy
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
// #include <core/scoring/methods/RG_Energy.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>

//Auto Headers
#include <core/conformation/Atom.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <ObjexxFCL/FArray.fwd.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class RG_Energy_Tests : public CxxTest::TestSuite {

	public:

	Pose pose;
	RG_Energy_Fast rg_energy_fast;
	core::chemical::ResidueTypeSetCAP rsd_set;
	core::scoring::ScoreFunctionOP scorefxn;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		//extern int command_line_argc; extern char ** command_line_argv;
		using namespace core;
		core_init();

		// scorefxn isn't really used, but necessary for call to finalize_total_energy
		scorefxn = new core::scoring::ScoreFunction;

		// rg_energy = new RG_Energy;
		// rg_energy_fast = new RG_Energy_Fast;
	}

	// Shared finalization goes here.
	void tearDown() {
		scorefxn = 0;
	}

	// --------------- Test Cases --------------- //
	void test_eval_energy() {

		float const TOLERATED_ERROR( 1e-3 );
		float const rg_correct(15.0709);

		float rg_fast;

		EnergyMap emap;
		emap.zero();
		pose = create_test_in_pdb_pose();
		//core::import_pose::pose_from_pdb( pose, "core/io/test_in.pdb" );
		rg_energy_fast.finalize_total_energy( pose, *scorefxn, emap );
		rg_fast = emap[ rg ];
		TS_ASSERT_DELTA( rg_fast, rg_correct, TOLERATED_ERROR );
	} // test_eval_energy
};


