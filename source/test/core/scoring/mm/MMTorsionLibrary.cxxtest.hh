// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMTorsionLibrary.cxxtest.hh
/// @brief  test suite for core::scoring::mm::MMTorsionLibrary
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/mm/MMTorsionLibrary.hh>

// Project headers
#include <test/core/init_util.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Numeric headers
#include <numeric/conversions.hh>

// C++ headers
#include <string>

//Auto Headers
#include <core/chemical/MMAtomTypeSet.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.hh>


using namespace core;
using namespace core::scoring::mm;
using namespace core::chemical;

// --------------- Test Class --------------- //

class MMTorsionLibraryTests : public CxxTest::TestSuite {

public:

	MMTorsionLibraryOP mmtorsionlibrary;
	MMAtomTypeSetOP mmatomtypeset;
	mm_torsion_param_set* params1;
	mm_torsion_param_set* params2;
	mm_torsion_param_set* params3;

	// --------------- Suite-level Fixture --------------- //

	MMTorsionLibraryTests() {
		core_init();

		// Only want to read in copy of the library once not for each test so init here in ctor

		// init the mmatomtypeset
		mmatomtypeset = MMAtomTypeSetOP( new MMAtomTypeSet );
		mmatomtypeset->read_file( "core/chemical/mm_atom_properties.txt" );

		// init the mmtorsionlibrary
		mmtorsionlibrary = MMTorsionLibraryOP( new MMTorsionLibrary( "core/scoring/mm/mm_torsion_params.txt" , mmatomtypeset ) );
	}

	virtual ~MMTorsionLibraryTests() {}

	static MMTorsionLibraryTests* createSuite() {
		return new MMTorsionLibraryTests();
	}

	static void destroySuite( MMTorsionLibraryTests *suite ) {
		delete suite;
	}

	// --------------- Fixtures --------------- //

	void setUp() {
		params1 = new mm_torsion_param_set( 0.2000, 1, numeric::conversions::radians( 180.00 ) );
		params2 = new mm_torsion_param_set( 0.2000, 2, numeric::conversions::radians( 180.00 ) );
		params3 = new mm_torsion_param_set( 0.2000, 3, numeric::conversions::radians( 180.00 ) );
	}

	void tearDown() {
		delete params1;
		delete params2;
		delete params3;
	}

	// ------------- Helper Function ------------- //

	mm_torsion_atom_quad make_quad( std::string s1,  std::string s2, std::string s3, std::string s4 ) {
		return mm_torsion_atom_quad( mmatomtypeset->atom_type_index( s1 ),
			mmatomtypeset->atom_type_index( s2 ),
			mmatomtypeset->atom_type_index( s3 ),
			mmatomtypeset->atom_type_index( s4 ) );
	}

	// --------------- Test Cases --------------- //

	void test_lookup() {

		mm_torsion_library_citer_pair mtlcp;
		mm_torsion_library_citer mtlc;

		// test lookup fully assigned by string
		// forward
		mm_torsion_atom_quad A = make_quad( "C", "CA", "CC", "CP1" );
		mtlcp = mmtorsionlibrary->lookup( "C", "CA", "CC", "CP1" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, A ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, A );
		// backward
		mtlcp =  mmtorsionlibrary->lookup( "CP1", "CC", "CA", "C" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, A ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, A );


		// test lookup fully assigned by string multiparam
		// forward
		mm_torsion_atom_quad B = make_quad( "CPH1", "CP3", "CC", "C" );
		mtlcp =  mmtorsionlibrary->lookup( "CPH1", "CP3", "CC", "C" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, B ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, B ); TS_ASSERT_EQUALS( mtlc->second, *params2 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, B ); TS_ASSERT_EQUALS( mtlc->second, *params3 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, B );
		// backward
		mtlcp =  mmtorsionlibrary->lookup( "C", "CC", "CP3", "CPH1" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, B ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, B ); TS_ASSERT_EQUALS( mtlc->second, *params2 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, B ); TS_ASSERT_EQUALS( mtlc->second, *params3 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, B );


		// test lookup explicite wildcard by string
		// forward
		mm_torsion_atom_quad C = make_quad( "X", "CP2", "CP3", "X" );
		mtlcp =  mmtorsionlibrary->lookup( "X", "CP2", "CP3", "X" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, C ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, C );
		// backward
		mtlcp =  mmtorsionlibrary->lookup( "X", "CP3", "CP2", "X" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, C ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, C );

		// test lookup explicite wildcard by string multiparam
		// forward
		mm_torsion_atom_quad D = make_quad( "X", "C", "CA", "X" );
		mtlcp =  mmtorsionlibrary->lookup( "X", "C", "CA", "X" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, D ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, D ); TS_ASSERT_EQUALS( mtlc->second, *params2 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, D ); TS_ASSERT_EQUALS( mtlc->second, *params3 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, D );
		// backward
		mtlcp =  mmtorsionlibrary->lookup( "X", "CA", "C", "X" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, D ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, D ); TS_ASSERT_EQUALS( mtlc->second, *params2 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, D ); TS_ASSERT_EQUALS( mtlc->second, *params3 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, D );

		// test lookup wildcard by string
		// forward
		mm_torsion_atom_quad E = make_quad( "X", "CP2", "CP3", "X" );
		mtlcp =  mmtorsionlibrary->lookup( "CC", "CP2", "CP3", "CC" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, E ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, E );
		// backward
		mtlcp =  mmtorsionlibrary->lookup( "CC", "CP3", "CP2", "CC" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, E ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, E );

		// test lookup wildcard by string multiparam
		// forward
		mm_torsion_atom_quad F = make_quad( "X", "C", "CA", "X" );
		mtlcp =  mmtorsionlibrary->lookup( "CC", "C", "CA", "CC" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, F ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, F ); TS_ASSERT_EQUALS( mtlc->second, *params2 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, F ); TS_ASSERT_EQUALS( mtlc->second, *params3 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, F );
		// backward
		mtlcp =  mmtorsionlibrary->lookup( "CC", "CA", "C", "CC" );
		mtlc = mtlcp.first;
		TS_ASSERT_EQUALS(  mtlc->first, F ); TS_ASSERT_EQUALS( mtlc->second, *params1 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, F ); TS_ASSERT_EQUALS( mtlc->second, *params2 ); mtlc++;
		TS_ASSERT_EQUALS(  mtlc->first, F ); TS_ASSERT_EQUALS( mtlc->second, *params3 ); mtlc++;
		TS_ASSERT_DIFFERS( mtlc->first, F );
	}
};
