// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/AtomTypeSet.cxxtest.hh
/// @brief  test suite for core::chemical::AtomTypeSet.cc
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/ElementSet.hh>
#include <core/chemical/Element.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

using namespace core;
using namespace core::chemical;

// --------------- Test Class --------------- //

class ElementSetTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	ElementSetOP element_set;
	double delta_percent;

	// --------------- Suite-level Fixture --------------- //

	AtomTypeSetTests() {
		//std::string commandline = "core.test -mute all";
		//initialize_from_commandline( commandline );
		core_init();


		// Want to read the properties file in only once for all the tests in this suite
		// so do it in the constructor for the test suite.

		// note this reads the element properties in the unit test directory
		// not the rosetta_database
		element_set = new ElementSet;
		element_set->read_file( "core/chemical/element_properties.txt" );
	}

	virtual ~ElementSetTests() {}

	static ElementSetTests *createSuite() {
		return new ElementSetTests();
	}

	static void destroySuite( ElementSetTests *suite ) {
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

	void test_n_elements() {
		TS_ASSERT_EQUALS( element_set->n_elements(), 17 );
	}

	void test_atom_type_index() {
		TS_ASSERT_EQUALS( element_set->element_index("H"), 1 );
		TS_ASSERT_EQUALS( element_set->element_index("C"), 2 );
		TS_ASSERT_EQUALS( element_set->element_index("N"), 3 );
		TS_ASSERT_EQUALS( element_set->element_index("O"), 4 );
		TS_ASSERT_EQUALS( element_set->element_index("F"), 5 );
		TS_ASSERT_EQUALS( element_set->element_index("NA"), 6 );
		TS_ASSERT_EQUALS( element_set->element_index("MG"), 7 );
	}

	void test_contains_element_type() {
		TS_ASSERT( element_set->contains_element_type("P"))
		TS_ASSERT( element_set->contains_element_type("S"))
		TS_ASSERT( element_set->contains_element_type("CL"))
		TS_ASSERT( !element_set->contains_element_type("FAKE"))
	}
	void test_operator_brackets() {
		// Test the parsing of the following line from the properties file
		//NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME
		//OH       O    1.5500    0.1591   -6.7700    3.5000   10.8000 ACCEPTOR SP3_HYBRID DONOR
		// This test is also testing some of the AtomType methods!
		ElementSet & e_set = *element_set; // just for convenience
		TS_ASSERT_EQUALS( element_set[1].symbol(), "H");
		TS_ASSERT_EQUALS( element_set[1].name(), "Hydrogen");
		TS_ASSERT_DELTA( element_set[1].weight(), 1.008, delta_percent);
		TS_ASSERT_DELTA( element_set[1].mass(), 1.008, delta_percent);
	}

};


