// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/nv/NVlookup.cxxtest.hh
/// @brief  test suite for core::scoring::NVlookup.cc
/// @author Sam DeLuca

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/nv/NVlookup.hh>

#include <platform/types.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/AA.hh>


//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::scoring;
using namespace core::scoring::nv;

class NVlookupTests : public CxxTest::TestSuite {

public:

	NVlookup const * nv_lookup;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		using namespace std;

		using namespace core;
		core_init();

		nv_lookup = & ScoringManager::get_instance()->get_NVLookupTable();

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_lookup_values() {

		float const TOLERATED_ERROR = 0.0001;


		EnergyMap emap;

		Real W_arbitrary_correct = 0.324105;
		Real W_zero_correct = -1.57599;
		Real W_one_correct = 5.19647;
		Real A_arbitrary_correct = 0.144313;

		//look up a few values out of the emap to make sure it's loading things properly
		Real W_arbitrary = nv_lookup->get_potentials(core::chemical::aa_trp,0.43);
		Real W_zero = nv_lookup->get_potentials(core::chemical::aa_trp,0.0);
		Real W_one = nv_lookup->get_potentials(core::chemical::aa_trp,1.0);
		Real A_arbitrary = nv_lookup->get_potentials(core::chemical::aa_ala,0.43);

		TS_ASSERT_DELTA(W_arbitrary,W_arbitrary_correct , TOLERATED_ERROR );
		TS_ASSERT_DELTA(W_zero,W_zero_correct,TOLERATED_ERROR);
		TS_ASSERT_DELTA(W_one,W_one_correct,TOLERATED_ERROR);
		TS_ASSERT_DELTA(A_arbitrary,A_arbitrary_correct,TOLERATED_ERROR);
	}

};


