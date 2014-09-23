// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleLibrary.cxxtest.hh
/// @brief  test suite for core::scoring::mm::MMBondAngleLibrary
/// @author Colin A. Smith (colin.smith@ucsf.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
// AUTO-REMOVED #include <basic/database/open.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>
// AUTO-REMOVED #include <core/types.hh>

// Project headers
#include <test/core/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Numeric headers
#include <numeric/conversions.hh>

// C++ headers
// AUTO-REMOVED #include <iostream>
#include <string>

//Auto Headers
#include <core/chemical/MMAtomTypeSet.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.hh>


using namespace core;
using namespace core::scoring::mm;
using namespace core::chemical;

// --------------- Test Class --------------- //

class MMBondAngleLibraryTests : public CxxTest::TestSuite {

public:

	MMBondAngleLibraryOP mmbondanglelibrary;
	MMAtomTypeSetOP mmatomtypeset;
	mm_bondangle_param_set* params1;
	mm_bondangle_param_set* params2;
	mm_bondangle_param_set* params3;

	// --------------- Suite-level Fixture --------------- //

	MMBondAngleLibraryTests() {
		core_init();

		// Only want to read in copy of the library once not for each test so init here in ctor

		// init the mmatomtypeset
		mmatomtypeset = MMAtomTypeSetOP( new MMAtomTypeSet );
		mmatomtypeset->read_file( "core/chemical/mm_atom_properties.txt" );

		// init the mmbondanglelibrary
		mmbondanglelibrary = MMBondAngleLibraryOP( new MMBondAngleLibrary( "core/scoring/mm/par_all27_prot_na.prm" , MMAtomTypeSetAP( mmatomtypeset ) ) );
	}

	virtual ~MMBondAngleLibraryTests() {}

	static MMBondAngleLibraryTests* createSuite() {
		return new MMBondAngleLibraryTests();
	}

	static void destroySuite( MMBondAngleLibraryTests *suite ) {
		delete suite;
	}

	// --------------- Fixtures --------------- //

	void setUp() {
		params1 = new mm_bondangle_param_set( 0.2000, numeric::conversions::radians( 180.00 ) );
		params2 = new mm_bondangle_param_set( 0.4000, numeric::conversions::radians( 180.00 ) );
		params3 = new mm_bondangle_param_set( 0.6000, numeric::conversions::radians( 180.00 ) );
	}

	void tearDown() {
		delete params1;
		delete params2;
		delete params3;
	}

	// ------------- Helper Function ------------- //

	mm_bondangle_atom_tri make_tri( std::string s1,  std::string s2, std::string s3 ) {
		return mm_bondangle_atom_tri( mmatomtypeset->atom_type_index( s1 ),
			mmatomtypeset->atom_type_index( s2 ),
			mmatomtypeset->atom_type_index( s3 ) );
	}

	// --------------- Test Cases --------------- //

		void test_lookup() {

		mm_bondangle_library_citer_pair mbalcp;
		mm_bondangle_library_citer mbalc;

		// test lookup fully assigned by string
		// forward
		mm_bondangle_atom_tri A = make_tri( "C", "CA", "CP1" );
		mbalcp = mmbondanglelibrary->lookup( "C", "CA", "CP1" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, A ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, A );
		// backward
		mbalcp =  mmbondanglelibrary->lookup( "CP1", "CA", "C" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, A ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, A );


		// test lookup fully assigned by string multiparam
		// forward
		mm_bondangle_atom_tri B = make_tri( "CPH1", "CP3", "C" );
		mbalcp =  mmbondanglelibrary->lookup( "CPH1", "CP3", "C" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, B ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, B ); TS_ASSERT_EQUALS( mbalc->second, *params2 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, B ); TS_ASSERT_EQUALS( mbalc->second, *params3 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, B );
		// backward
		mbalcp =  mmbondanglelibrary->lookup( "C", "CP3", "CPH1" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, B ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, B ); TS_ASSERT_EQUALS( mbalc->second, *params2 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, B ); TS_ASSERT_EQUALS( mbalc->second, *params3 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, B );


		// test lookup explicite wildcard by string
		// forward
		mm_bondangle_atom_tri C = make_tri( "X", "CP2", "X" );
		mbalcp =  mmbondanglelibrary->lookup( "X", "CP2", "X" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, C ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, C );
		// backward
		mbalcp =  mmbondanglelibrary->lookup( "X", "CP2", "X" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, C ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, C );

		// test lookup explicite wildcard by string multiparam
		// forward
		mm_bondangle_atom_tri D = make_tri( "X", "C", "X" );
		mbalcp =  mmbondanglelibrary->lookup( "X", "C", "X" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, D ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, D ); TS_ASSERT_EQUALS( mbalc->second, *params2 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, D ); TS_ASSERT_EQUALS( mbalc->second, *params3 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, D );
		// backward
		mbalcp =  mmbondanglelibrary->lookup( "X", "C", "X" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, D ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, D ); TS_ASSERT_EQUALS( mbalc->second, *params2 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, D ); TS_ASSERT_EQUALS( mbalc->second, *params3 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, D );

		// test lookup wildcard by string
		// forward
		mm_bondangle_atom_tri E = make_tri( "X", "CP2", "X" );
		mbalcp =  mmbondanglelibrary->lookup( "CC", "CP2", "CA" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, E ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, E );
		// backward
		mbalcp =  mmbondanglelibrary->lookup( "CA", "CP2", "CC" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, E ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, E );

		// test lookup wildcard by string multiparam
		// forward
		mm_bondangle_atom_tri F = make_tri( "X", "C", "X" );
		mbalcp =  mmbondanglelibrary->lookup( "CC", "C", "CA" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, F ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, F ); TS_ASSERT_EQUALS( mbalc->second, *params2 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, F ); TS_ASSERT_EQUALS( mbalc->second, *params3 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, F );
		// backward
		mbalcp =  mmbondanglelibrary->lookup( "CA", "C", "CC" );
		mbalc = mbalcp.first;
		TS_ASSERT_EQUALS(  mbalc->first, F ); TS_ASSERT_EQUALS( mbalc->second, *params1 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, F ); TS_ASSERT_EQUALS( mbalc->second, *params2 ); mbalc++;
		TS_ASSERT_EQUALS(  mbalc->first, F ); TS_ASSERT_EQUALS( mbalc->second, *params3 ); mbalc++;
		TS_ASSERT_DIFFERS( mbalc->first, F );
	}
};
