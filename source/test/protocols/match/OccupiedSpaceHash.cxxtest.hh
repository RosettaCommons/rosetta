// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/OccupiedSpaceHash.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>
#include <protocols/match/OccupiedSpaceHash.hh>

// C++ headers
#include <string>
#include <iostream>

#include <numeric/HomogeneousTransform.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map_fwd.hpp>


// --------------- Test Class --------------- //

class OccupiedSpaceHashTests : public CxxTest::TestSuite {

	public:

		typedef core::Vector Vector;
		typedef core::Size   Size;
		typedef core::Real   Real;
		typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
		typedef numeric::HomogeneousTransform< Real >    HTReal;
	  typedef numeric::geometry::hashing::Real6 Real6;

	private:
		Vector lower;
		Vector upper;
		BoundingBox bb;

		numeric::geometry::hashing::Real6 pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ;

	public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		lower = Vector( 12.5, 16.25, 4.25 );
		upper = Vector( 15.5, 20,    8.5  );
		bb = BoundingBox( lower, upper );

		pA[ 1 ] = 13.6; pA[ 2 ] = 19.4; pA[ 3 ] = 5.3; pA[ 4 ] = 50;   pA[ 5 ] = 123;  pA[ 6 ] = 76;
		pB[ 1 ] = 13.6; pB[ 2 ] = 19.3; pB[ 3 ] = 5.3; pB[ 4 ] = 54;   pB[ 5 ] = 123;  pB[ 6 ] = 78;
		pC[ 1 ] = 13.4; pC[ 2 ] = 19.4; pC[ 3 ] = 5.3; pC[ 4 ] = 53;   pC[ 5 ] = 133;  pC[ 6 ] = 76;
		pD[ 1 ] = 13.6; pD[ 2 ] = 19.4; pD[ 3 ] = 5.3; pD[ 4 ] = 70;   pD[ 5 ] = 123;  pD[ 6 ] = 76;
		pE[ 1 ] = 13.6; pE[ 2 ] = 19.4; pE[ 3 ] = 5.3; pE[ 4 ] = 75;   pE[ 5 ] = 120;  pE[ 6 ] = 76;
		pF[ 1 ] = 13.7; pF[ 2 ] = 19.4; pF[ 3 ] = 5.4; pF[ 4 ] = 75;   pF[ 5 ] = 123;  pF[ 6 ] = 76;
		pG[ 1 ] = 13.8; pG[ 2 ] = 19.4; pG[ 3 ] = 5.4; pG[ 4 ] = 75;   pG[ 5 ] = 123;  pG[ 6 ] = 73;
		pH[ 1 ] = 13.9; pH[ 2 ] = 19.4; pH[ 3 ] = 5.4; pH[ 4 ] = 75;   pH[ 5 ] = 123;  pH[ 6 ] = 73;
		pI[ 1 ] = 14.0; pI[ 2 ] = 19.4; pI[ 3 ] = 5.5; pI[ 4 ] = 75;   pI[ 5 ] = 123;  pI[ 6 ] = 76;
		pJ[ 1 ] = 16.0; pJ[ 2 ] = 19.4; pJ[ 3 ] = 5.5; pJ[ 4 ] = 75;   pJ[ 5 ] = 123;  pJ[ 6 ] = 76;
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //
	void test_occupied_space_hash_ctor() {
		using namespace protocols::match;

		OccupiedSpaceHash space;
		space.set_bounding_box( bb );
		space.set_uniform_xyz_bin_width( 0.25 );
		space.set_uniform_euler_angle_bin_width( 10 );

		space.initialize();

		/// After it's construction, the OccSpaceHash knows that it's merely performing
		/// bounding-box checks for the points its given.
		TS_ASSERT( space.match_possible_for_hit_geometry( pA ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pB ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pC ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pD ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pE ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pF ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pG ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pH ) );
		TS_ASSERT( space.match_possible_for_hit_geometry( pI ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pJ ) );

	}

	void test_occupied_space_hash_note_hit_geometry() {
		using namespace protocols::match;
		OccupiedSpaceHash space;
		space.set_bounding_box( bb );
		space.set_uniform_xyz_bin_width( 0.25 );
		space.set_uniform_euler_angle_bin_width( 10 );

		space.initialize();

		space.insert_hit_geometry( pA );

		TS_ASSERT(   space.match_possible_for_hit_geometry( pA ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pB ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pC ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pD ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pE ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pF ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pG ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pH ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pI ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pJ ) );

		space.insert_hit_geometry( pF );

		TS_ASSERT(   space.match_possible_for_hit_geometry( pA ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pB ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pC ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pD ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pE ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pF ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pG ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pH ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pI ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pJ ) );

	}

	void test_occupied_space_hash_note_hit_geometry2() {
		using namespace protocols::match;
		OccupiedSpaceHash space;
		space.set_bounding_box( bb );
		space.set_uniform_xyz_bin_width( 0.25 );
		space.set_uniform_euler_angle_bin_width( 10 );

		space.initialize();

		space.insert_hit_geometry( pB );

		TS_ASSERT(   space.match_possible_for_hit_geometry( pA ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pB ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pC ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pD ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pE ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pF ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pG ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pH ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pI ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pJ ) );

		space.insert_hit_geometry( pI );

		TS_ASSERT(   space.match_possible_for_hit_geometry( pA ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pB ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pC ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pD ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pE ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pF ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pG ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pH ) );
		TS_ASSERT(   space.match_possible_for_hit_geometry( pI ) );
		TS_ASSERT( ! space.match_possible_for_hit_geometry( pJ ) );

	}

};
