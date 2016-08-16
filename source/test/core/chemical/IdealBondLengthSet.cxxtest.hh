// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/AtomTypeSet.cxxtest.hh
/// @brief  test suite for core::chemical::AtomTypeSet.cc
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/IdealBondLengthSet.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

using namespace core;
using namespace core::chemical;

// --------------- Test Class --------------- //

class IdealBondLengthSetTests : public CxxTest::TestSuite {

	public:

	// Shared data ideal_bond_lengths go here.
	IdealBondLengthSetOP ideal_bond_length_set;
	double delta_percent;

	// --------------- Suite-level Fixture --------------- //

	IdealBondLengthSetTests() {
		//std::string commandline = "core.test -mute all";
		//initialize_from_commandline( commandline );
		core_init();


		// Want to read the properties file in only once for all the tests in this suite
		// so do it in the constructor for the test suite.

		// note this reads the ideal_bond_length properties in the unit test directory
		// not the rosetta_database
		ideal_bond_length_set = new IdealBondLengthSet;
		ideal_bond_length_set->read_file( "core/chemical/ideal_bond_lengths.txt" );
	}

	virtual ~IdealBondLengthSetTests() {}

	static IdealBondLengthSetTests *createSuite() {
		return new IdealBondLengthSetTests();
	}

	static void destroySuite( IdealBondLengthSetTests *suite ) {
		delete suite;
	}


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. The fixture above
	// gets constructed once before all the tests in this test suite are run.

	// Shared initialization goes here.
	void setUp() {
		delta_percent = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {	}


	// --------------- Test Cases --------------- //

	void test_contains_bond_length() {
		TS_ASSERT( ideal_bond_length_set->contains_bond_length("aroC", "CAbb"));
		TS_ASSERT( ideal_bond_length_set->contains_bond_length("CAbb", "aroC"));
		TS_ASSERT( ideal_bond_length_set->contains_bond_length("VIRT", "VIRT"));

		TS_ASSERT(! ideal_bond_length_set->contains_bond_length("Fake", "Made-up"));
	}

	void test_get_bond_length() {
		TS_ASSERT_DELTA( ideal_bond_length_set->get_bond_length("aroC", "CAbb"), 1.529832, delta_percent);
		TS_ASSERT_DELTA( ideal_bond_length_set->get_bond_length("CAbb", "aroC"), 1.529832, delta_percent);
		TS_ASSERT_DELTA( ideal_bond_length_set->get_bond_length("VIRT", "VIRT"), 1.569086, delta_percent);

		//TS_ASSERT_THROWS(ideal_bond_length_set->get_bond_length("Fake", "Made-up");// code exits, doesn' throw exception uet
		
		TS_ASSERT_EQUALS(
				ideal_bond_length_set->get_bond_length("CH2", "Narg")
				ideal_bond_length_set->get_bond_length("Narg", "CH2")
		)
	}

};


