// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ContactOrderEnergy.cxxtest.hh
/// @brief  test suite for ContactOrderEnergy
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/ContactOrderEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <core/chemical/ChemicalManager.hh>
#include <test/core/init_util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/silent/SilentFileData.hh>

//Auto Headers
#include <core/io/silent/EnergyNames.fwd.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class ContactOrderEnergy_Tests : public CxxTest::TestSuite {

public:

	Pose pose;
	ContactOrderEnergy co_energy;
	core::io::silent::SilentFileData sfd;
	core::chemical::ResidueTypeSetCOP rsd_set;
	core::scoring::ScoreFunctionOP scorefxn;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		using namespace std;

		//extern int command_line_argc; extern char ** command_line_argv;
		using namespace core;
		core_init_with_additional_options( "-in::file::silent_struct_type protein" );

		// correct answers taken from rosetta++ v19429
		sfd.read_file( "core/scoring/methods/score3_in.silent_out" );
		// contact order calculations in rosetta++ used centroids
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		// scorefxn isn't really used, but necessary for call to finalize_total_energy
		scorefxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	}

	// Shared finalization goes here.
	void tearDown() {
		scorefxn.reset();
	}

	// --------------- Test Cases --------------- //
	void test_eval_energy() {

		// TOLERATED_ERROR is set to the largest deviation from calculating contact order in mini and r++
		// as of 8/13/08. While this allowed deviation is large, out of the ~1k structures in
		// score3_in.silent_out there are only 17 with an error greater than 0.1.
		float const TOLERATED_ERROR = 0.317;

		float co, silent_co;
		EnergyMap emap;
		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
			emap.zero();
			iter->fill_pose( pose, *rsd_set );
			silent_co = iter->get_energy( "co" );

			co_energy.finalize_total_energy( pose, *scorefxn, emap );
			co = emap[ core::scoring::co ];
			TS_ASSERT_DELTA( co, silent_co, TOLERATED_ERROR );
		}
	} // test_eval_energy()
}; // ContactOrderEnergy_Tests


