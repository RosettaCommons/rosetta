// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/MMAtomTypeSet.cxxtest.hh
/// @brief  test suite for core::chemical::MMAtomTypeSet
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers, for debugging your tests

//Auto Headers


using namespace core;
using namespace core::chemical;

// --------------- Test Class --------------- //

class MMAtomTypeSetTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	MMAtomTypeSetOP mmatomtypeset;

	// --------------- Suite-level Fixture --------------- //

	MMAtomTypeSetTests() {
		//std::string commandline = "core.test -mute all";
		//initialize_from_commandline( commandline );
		core_init();

		// Want to read the properties file in only once for all the tests in this suite
		// so do it in the constructor for the test suite.
		mmatomtypeset = MMAtomTypeSetOP( new MMAtomTypeSet() );
		mmatomtypeset->read_file( "core/chemical/mm_atom_properties.txt" );
	}

	virtual ~MMAtomTypeSetTests() {}

	static MMAtomTypeSetTests* createSuite() {
		return new MMAtomTypeSetTests();
	}

	static void destroySuite( MMAtomTypeSetTests *suite ) {
		delete suite;
	}

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	void test_n_atomtypes() {
		TS_ASSERT_EQUALS( mmatomtypeset->n_atomtypes(), 9 );
	}

	void test_contains_atom_type() {
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("C"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("CA"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("CC"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("CP1"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("CP2"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("CP3"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("CPH1"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("X"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("VIRT"), true);
		TS_ASSERT_EQUALS(mmatomtypeset->contains_atom_type("N"), false);
	}

	void test_atom_type_index() {
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("C"),    1 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("CA"),   2 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("CC"),   3 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("CP1"),  4 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("CP2"),  5 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("CP3"),  6 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("CPH1"), 7 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("X"),    8 );
		TS_ASSERT_EQUALS( mmatomtypeset->atom_type_index("VIRT"), 9 );
	}

	void test_operator_brackets() {
		// This tests the MMAtomType as well
		MMAtomType C("C");
		TS_ASSERT_EQUALS( (*mmatomtypeset)[1].name(), C.name() );
	}


};
