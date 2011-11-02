// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/orbitals/OrbitalTypeSet.cxxtest.hh
/// @brief  test suite for core::chemical::orbitals::OrbitalTypeSet.cc
/// @author Steven Combs


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh> // for some reason, AtomTypeSet.fwd.hh does not #include owning_ptr.hh

// C++ headers, for debugging your tests
// AUTO-REMOVED #include <iostream>

//Auto Headers
#include <utility/vector1.hh>


using namespace core;
using namespace core::chemical;

// --------------- Test Class --------------- //

class OrbitalTypeSetTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	orbitals::OrbitalTypeSetOP orbitaltypeset;
	double delta_percent;

	// --------------- Suite-level Fixture --------------- //

	OrbitalTypeSetTests() {
		//std::string commandline = "core.test -mute all";
		//initialize_from_commandline( commandline );
		core_init();


		// Want to read the properties file in only once for all the tests in this suite
		// so do it in the constructor for the test suite.
		orbitaltypeset = new orbitals::OrbitalTypeSet( "core/chemical/orbitals/");
		//atomtypeset->read_file( "core/chemical/atom_properties.txt" );
	}

	virtual ~OrbitalTypeSetTests() {}

	static OrbitalTypeSetTests *createSuite() {
		return new OrbitalTypeSetTests();
	}

	static void destroySuite( OrbitalTypeSetTests *suite ) {
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


	void test_orbital_type_index() {
		TS_ASSERT_EQUALS( orbitaltypeset->orbital_type_index("C.pi.sp2"), 1 );
		TS_ASSERT_EQUALS( orbitaltypeset->orbital_type_index("N.pi.sp2"), 2 );
		TS_ASSERT_EQUALS( orbitaltypeset->orbital_type_index("N.p.sp2"), 3 );
		TS_ASSERT_EQUALS( orbitaltypeset->orbital_type_index("O.pi.sp2"), 4 );
		TS_ASSERT_EQUALS( orbitaltypeset->orbital_type_index("O.p.sp2"), 5 );
		TS_ASSERT_EQUALS( orbitaltypeset->orbital_type_index("O.p.sp3"), 6 );
		TS_ASSERT_EQUALS( orbitaltypeset->orbital_type_index("S.p.sp3"), 7 );
	}

	void test_operator_brackets() {
		// Test the parsing of the following line from the properties file
		// This tests the OrbitalType functions! This means I do not have to write orbitaltype functions!!!!
		orbitals::OrbitalTypeSet & orbitalsetr = *orbitaltypeset;
		TS_ASSERT_DELTA( orbitalsetr[4].distance(), 1.100, delta_percent);
		TS_ASSERT_EQUALS( orbitalsetr[4].hybridization(), "sp2");
		TS_ASSERT_EQUALS( orbitalsetr[4].orbital_name(), "pi");
		TS_ASSERT_EQUALS( orbitalsetr[4].name(),"O.pi.sp2");
		TS_ASSERT_EQUALS( orbitalsetr[6].hybridization(), "sp3");
		TS_ASSERT_EQUALS( orbitalsetr[6].orbital_name(), "p" );
		TS_ASSERT_EQUALS( orbitalsetr[6].name(), "O.p.sp3");
	}


};



