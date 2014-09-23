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
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

using namespace core;
using namespace core::chemical;

// --------------- Test Class --------------- //

class AtomTypeSetTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	AtomTypeSetOP atomtypeset;
	double delta_percent;


	// --------------- Suite-level Fixture --------------- //

	AtomTypeSetTests() {
		//std::string commandline = "core.test -mute all";
		//initialize_from_commandline( commandline );
		core_init();


		// Want to read the properties file in only once for all the tests in this suite
		// so do it in the constructor for the test suite.

		// note this reads the atom properties in the unit test directory
		// not the rosetta_database
		atomtypeset = AtomTypeSetOP( new AtomTypeSet( "core/chemical/") );
		//atomtypeset->read_file( "core/chemical/atom_properties.txt" );
	}

	virtual ~AtomTypeSetTests() {}

	static AtomTypeSetTests *createSuite() {
		return new AtomTypeSetTests();
	}

	static void destroySuite( AtomTypeSetTests *suite ) {
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
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	void test_atom_type_set_name() {
		TS_ASSERT_EQUALS(atomtypeset->name(), "chemical");
	}

	void test_n_atomtypes() {
		TS_ASSERT_EQUALS( atomtypeset->n_atomtypes(), 7 );
	}

	void test_atom_type_index() {
		TS_ASSERT_EQUALS( atomtypeset->atom_type_index("CNH2"), 1 );
		TS_ASSERT_EQUALS( atomtypeset->atom_type_index("Ntrp"), 2 );
		TS_ASSERT_EQUALS( atomtypeset->atom_type_index("Nhis"), 3 );
		TS_ASSERT_EQUALS( atomtypeset->atom_type_index("OH"), 4 );
		TS_ASSERT_EQUALS( atomtypeset->atom_type_index("ONH2"), 5 );
		TS_ASSERT_EQUALS( atomtypeset->atom_type_index("OOC"), 6 );
		TS_ASSERT_EQUALS( atomtypeset->atom_type_index("Hpol"), 7 );
	}

	void test_operator_brackets() {
		// Test the parsing of the following line from the properties file
		//NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME
		//OH       O    1.5500    0.1591   -6.7700    3.5000   10.8000 ACCEPTOR SP3_HYBRID DONOR
		// This test is also testing some of the AtomType methods!
		AtomTypeSet & atomsetr = *atomtypeset;
		TS_ASSERT_DELTA( atomsetr[4].lk_lambda(), 3.5000, delta_percent);
		TS_ASSERT_DELTA( atomsetr[4].lk_dgfree(), -6.7700, delta_percent);
		TS_ASSERT_DELTA( atomsetr[4].lk_volume(), 10.8000, delta_percent);
		TS_ASSERT_DELTA( atomsetr[4].lj_radius(), 1.5500, delta_percent);
		TS_ASSERT_DELTA( atomsetr[4].lj_wdepth(), 0.1591, delta_percent);
		TS_ASSERT( atomsetr[4].is_acceptor() );
		TS_ASSERT( atomsetr[4].is_donor() );
		TS_ASSERT( !(atomsetr[4].is_hydrogen()) );
		TS_ASSERT( atomsetr[4].is_heavyatom() );
		TS_ASSERT( !(atomsetr[4].is_h2o()) );
		TS_ASSERT_EQUALS( atomsetr[4].hybridization(), SP3_HYBRID );
	}

	// Should there be more tests for AtomType set_parameter() and set_property()??

};


