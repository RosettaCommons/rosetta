// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/gasteiger/GasteigerAtomTypeSet.cxxtest.hh
/// @brief  test suite for core::chemical::gasteiger::GasteigerAtomTypeSet.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/Element.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

using namespace core;
using namespace core::chemical;
using namespace core::chemical::gasteiger;

// --------------- Test Class --------------- //

class GasteigerAtomTypeSetTests : public CxxTest::TestSuite {

	public:

	GasteigerAtomTypeSetCOP atom_type_set_;

	double delta_percent;

	// --------------- Suite-level Fixture --------------- //

	GasteigerAtomTypeSetTests() {
	}

	virtual ~GasteigerAtomTypeSetTests() {}

	static GasteigerAtomTypeSetTests *createSuite() {
		return new GasteigerAtomTypeSetTests();
	}

	static void destroySuite( GasteigerAtomTypeSetTests *suite ) {
		delete suite;
	}

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. The fixture above
	// gets constructed once before all the tests in this test suite are run.

	// Shared initialization goes here.
	void setUp() {
		core_init();

		atom_type_set_ = ChemicalManager::get_instance()->gasteiger_atom_type_set();
		delta_percent = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	void test_atom_type_retreval() {
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "C_TeTeTeTe" )->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "C_TeTeTeTe" )->get_number_hybrid_orbitals(), 4 );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "N_TrTrTrPi2" )->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "N_TrTrTrPi2" )->get_number_hybrid_orbitals(), 3 );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "H_S" )->get_name(), "H_S" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "H_S" )->get_number_hybrid_orbitals(), 0 );
	}

	void test_element_association() {
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "C_TeTeTeTe" )->get_element_type()->get_chemical_name(), "Carbon" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "N_TrTrTrPi2" )->get_element_type()->get_chemical_name(), "Nitrogen" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "H_S" )->get_element_type()->get_chemical_name(), "Hydrogen" );
	}

};


