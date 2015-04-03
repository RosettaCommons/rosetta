// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/gasteiger/GasteigerAtomTypeSet.cxxtest.hh
/// @brief  test suite for bond lengths
/// @author Steven Combs


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/bond_support.hh>
#include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

using namespace core;
using namespace core::chemical;
using namespace core::chemical::gasteiger;

// --------------- Test Class --------------- //

class BondLengthsTests : public CxxTest::TestSuite {

	public:

	GasteigerAtomTypeSetCOP atom_type_set_;


	// --------------- Suite-level Fixture --------------- //

	BondLengthsTests() {
	}

	virtual ~BondLengthsTests() {}

	static BondLengthsTests *createSuite() {
		return new BondLengthsTests();
	}

	static void destroySuite( BondLengthsTests *suite ) {
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
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	void test_bond_lengthsl() {
		TS_ASSERT_EQUALS(core::chemical::create_bond_length( *(atom_type_set_->atom_type( "C_TeTeTeTe" )), *(atom_type_set_->atom_type( "C_TeTeTeTe" )), core::chemical::SingleBond), 1.52);
		TS_ASSERT_EQUALS(core::chemical::create_bond_length( *(atom_type_set_->atom_type( "N_TrTrTrPi2" )), *(atom_type_set_->atom_type( "N_TrTrTrPi2" )), core::chemical::SingleBond), 1.42);
		TS_ASSERT_EQUALS(core::chemical::create_bond_length( *(atom_type_set_->atom_type( "N_TrTrTrPi2" )), *(atom_type_set_->atom_type( "H_S" )), core::chemical::SingleBond), 1.02);
		TS_ASSERT_EQUALS(core::chemical::create_bond_length( *(atom_type_set_->atom_type( "C_TeTeTeTe" )), *(atom_type_set_->atom_type( "H_S" )), core::chemical::SingleBond), 1.07);
		//std::cout << atom_type_set_->atom_type( "C_TeTeTeTe" )->get_atom_type_property(core::chemical::gasteiger::GasteigerAtomTypeData::CovalentRadiusSingleBond) << std::endl;
		//std::cout << core::chemical::create_bond_length( *(atom_type_set_->atom_type( "C_TeTeTeTe" )), *(atom_type_set_->atom_type( "C_TeTeTeTe" )), core::chemical::SingleBond) << std::endl;
		//std::cout << core::chemical::create_bond_length( *(atom_type_set_->atom_type( "N_TrTrTrPi2" )), *(atom_type_set_->atom_type( "N_TrTrTrPi2" )), core::chemical::SingleBond) << std::endl;
		//std::cout << core::chemical::create_bond_length( *(atom_type_set_->atom_type( "N_TrTrTrPi2" )), *(atom_type_set_->atom_type( "H_S" )), core::chemical::SingleBond) << std::endl;
		//std::cout << core::chemical::create_bond_length( *(atom_type_set_->atom_type( "C_TeTeTeTe" )), *(atom_type_set_->atom_type( "H_S" )), core::chemical::SingleBond) << std::endl;

/*		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "C_TeTeTeTe" )->get_name(), "C_TeTeTeTe" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "C_TeTeTeTe" )->get_number_hybrid_orbitals(), 4 );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "N_TrTrTrPi2" )->get_name(), "N_TrTrTrPi2" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "N_TrTrTrPi2" )->get_number_hybrid_orbitals(), 3 );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "H_S" )->get_name(), "H_S" );
		TS_ASSERT_EQUALS( atom_type_set_->atom_type( "H_S" )->get_number_hybrid_orbitals(), 0 );*/
	}


};


