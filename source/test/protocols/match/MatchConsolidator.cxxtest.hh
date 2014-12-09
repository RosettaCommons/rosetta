// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/Hasher.cxxtest.hh
/// @brief  test suite for protocols::match::SixDHasher
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Utility headers
#include <utility/fixedsizearray1.hh>
// AUTO-REMOVED #include <utility/exit.hh>

/// Project headers
#include <core/types.hh>
#include <protocols/match/output/MatchConsolidator.hh>

// C++ headers
// AUTO-REMOVED #include <string>
// AUTO-REMOVED #include <iostream>

//Auto Headers
#include <protocols/match/Hit.hh>
#include <utility/vector0.hh>
#include <boost/unordered/unordered_map_fwd.hpp>


using namespace protocols::match;
using namespace protocols::match::output;

// --------------- Test Class --------------- //

class MatchConsolidatorTests : public CxxTest::TestSuite {

public:

	typedef core::Vector Vector;
	typedef core::Size   Size;
	typedef core::Real   Real;


	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

  //FIX THIS ASAP

	// --------------- Test Cases --------------- //
	void test_BestMatchesCollection_not_full() {
		BestMatchesCollection bmc( 3, false );
		match m1( 2 ), m2( 2 );
		m1[ 1 ].first()[ 1 ] = 1; m2[ 1 ].first()[ 1 ] = 2;
		bmc.add_match( m1, 5 );
		bmc.add_match( m2, 6 );
		TS_ASSERT( bmc.n_kept_matches() == 2 );
		TS_ASSERT( bmc.kept_match( 1 )[ 1 ].first()[ 1 ] == 1 );
		TS_ASSERT( bmc.kept_match( 2 )[ 1 ].first()[ 1 ] == 2 );
		
	}

	void test_BestMatchesCollection_full() {
		BestMatchesCollection bmc( 3, false );
		match m1( 2 ), m2( 2 ), m3( 2 ), m4( 2 );
		m1[ 1 ].first()[ 1 ] = 1; m2[ 1 ].first()[ 1 ] = 2;
		m3[ 1 ].first()[ 1 ] = 3; m4[ 1 ].first()[ 1 ] = 4;
		bmc.add_match( m1, 5 );
		bmc.add_match( m2, 6 );
		bmc.add_match( m3, 4 );
		bmc.add_match( m4, 2 );
		TS_ASSERT( bmc.n_kept_matches() == 3 );
		TS_ASSERT( bmc.kept_match( 1 )[ 1 ].first()[ 1 ] == 1 );
		TS_ASSERT( bmc.kept_match( 2 )[ 1 ].first()[ 1 ] == 4 );
		TS_ASSERT( bmc.kept_match( 3 )[ 1 ].first()[ 1 ] == 3 );

	}


};
