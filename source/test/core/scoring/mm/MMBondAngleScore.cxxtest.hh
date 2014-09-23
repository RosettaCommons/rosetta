// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleScore.cxxtest.hh
/// @brief  test suite for core::scoring::mm::MMBondAngleScore
/// @author Colin A. Smith (colin.smith@ucsf.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>

// Unit headers
#include <core/scoring/mm/MMBondAngleScore.hh>
#include <core/types.hh>

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

class MMBondAngleScoreTests : public CxxTest::TestSuite {

public:

	MMBondAngleLibraryOP mmbondanglelibrary;
	MMAtomTypeSetOP mmatomtypeset;
	MMBondAngleScoreOP mmbondanglescore;

	// --------------- Suite-level Fixture --------------- //

	MMBondAngleScoreTests() {
		core_init();

		// Only want to read in copy of the library once not for each test so init here in ctor

		// init the mmatomtypeset
		mmatomtypeset = MMAtomTypeSetOP( new MMAtomTypeSet );
		mmatomtypeset->read_file( "core/chemical/mm_atom_properties.txt" );

		// init the mmbondanglelibrary
		mmbondanglelibrary = MMBondAngleLibraryOP( new MMBondAngleLibrary( "core/scoring/mm/par_all27_prot_na.prm" , MMAtomTypeSetAP( mmatomtypeset ) ) );

		// init the mmbondanglescore
		mmbondanglescore = MMBondAngleScoreOP( new MMBondAngleScore( *mmbondanglelibrary ) );
	 }

	virtual ~MMBondAngleScoreTests() {}

	static MMBondAngleScoreTests* createSuite() {
		return new MMBondAngleScoreTests();
	}

	static void destroySuite( MMBondAngleScoreTests *suite ) {
		delete suite;
	}

	// --------------- Fixtures --------------- //

	void setUp() {
	}

	void tearDown() { }

	// ------------- Helper Function ------------- //

	mm_bondangle_atom_tri make_tri( std::string s1,  std::string s2, std::string s3 ) {
		return mm_bondangle_atom_tri( mmatomtypeset->atom_type_index( s1 ),
			mmatomtypeset->atom_type_index( s2 ),
			mmatomtypeset->atom_type_index( s3 ) );
	}

	// --------------- Test Cases --------------- //

	void test_score() {

		test::UTracer UT("core/scoring/mm/MMBondAngleScoreTests.u");
		//std::ofstream UT("core/scoring/mm/MMBondAngleScoreTests.u");

		// make set of angles in radians
		Real angles[180];
		for( int i = 0; i<180; ++i )
			{
				angles[i] = numeric::conversions::radians( static_cast<Real>( i ) );
			}

		UT << "Single Parameter Scores:" << "\n";
		mm_bondangle_atom_tri A = make_tri( "C", "CA", "CP1" );
		for( int i = 0; i<180; ++i )
			{
				Real const score( mmbondanglescore->score( A, angles[i] ) );
				UT << numeric::conversions::degrees( angles[i] ) << "\t" << score << "\n";
			}

		UT << "Multiple Parameter Scores:" << "\n";
		mm_bondangle_atom_tri B = make_tri( "CPH1", "CP3", "C" );
		for(  int i = 0; i<180; ++i )
			{
				Real const score( mmbondanglescore->score( B, angles[i] ) );
				UT << numeric::conversions::degrees( angles[i] ) << "\t" << score << "\n";
			}
	}
};
