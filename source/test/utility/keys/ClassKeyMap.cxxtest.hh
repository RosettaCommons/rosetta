// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/ClassKeyMap.cxxtest.hh
/// @brief  ClassKeyMap unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/keys/ClassKeyMap.hh>


// --- set up the individual test cases
class KeyMapTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests for int
	void test_ClassKeyMap_general() {

		typedef  utility::keys::ClassKeyMap< int, int, float >  Map;

		// Create an instance and check starting conditions
		Map m;
		m.assign( 1, 123 ).assign( 2, 123 ).assign( 3, 123 );
		TS_ASSERT( m.size() == 3 );

		// Assign another key/value pair
		m.assign( 7, 77 );
		TS_ASSERT( m[ 7 ] == 77 );

		// Add a value with operator()
		m( 8 ) = 88;
		TS_ASSERT( m[ 8 ] == 88 );

		// Add a key
		m.add( 12 );
		TS_ASSERT( m.has( 12 ) );
		TS_ASSERT( m[ 12 ] == 0 );
		TS_ASSERT( m.size() == 6 );

		// Add another key
		m.add( 22 );
		TS_ASSERT( m.has( 22 ) );
		TS_ASSERT( m[ 22 ] == 0 );
		TS_ASSERT( m.size() == 7 );

		// Shrink to right-sized
		m.shrink();
		TS_ASSERT( m.capacity() == m.size() );
	}


	/// @brief Swap test
	void test_ClassKeyMap_swap() {

		typedef  utility::keys::ClassKeyMap< int, int, double >  Map;

		// Create an instance and check starting conditions
		Map m;
		m.assign( 3, 33 ).assign( 6, 33 ).assign( 9, 33 );
		TS_ASSERT( m.size() == 3 );

		// Creat a second ClassKeyMap
		Map w;
		w.assign( 3, 66 ).assign( 6, 66 ).assign( 9, 66 ).assign( 12, 66 ).assign( 15, 66 ).assign( 18, 66 );
		TS_ASSERT( w.size() == 6 );

		// Swap the two ClassKeyMaps
		swap( m, w );
		TS_ASSERT( m.size() == 6 );
		TS_ASSERT( w.size() == 3 );
		TS_ASSERT( m[ 3 ] == 66 );
		TS_ASSERT( w[ 9 ] == 33 );

		// Swap the two ClassKeyMaps back
		w.swap( m );
		TS_ASSERT( w.size() == 6 );
		TS_ASSERT( m.size() == 3 );
		TS_ASSERT( w[ 3 ] == 66 );
		TS_ASSERT( m[ 9 ] == 33 );
	}

};


