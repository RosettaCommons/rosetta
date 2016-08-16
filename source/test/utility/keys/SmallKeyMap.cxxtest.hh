// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/SmallKeyMap.cxxtest.hh
/// @brief  SmallKeyMap unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/keys/SmallKeyMap.hh>


class SmallKeyMapTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests
	void test_SmallKeyMap_general() {

		typedef  utility::keys::SmallKeyMap< int, int >  Map;

		// build initial map with 3 entries
		Map m;
		m.assign( 1, 123 ).assign( 2, 123 ).assign( 3, 123 );
		TS_ASSERT( m.size() == 3 );

		// add 4th entry
		m.assign( 7, 77 );
		TS_ASSERT( m[ 7 ] == 77 );

		// add 5th entry
		m( 8 ) = 88;

		// add 6th entry
		m.add( 12 );
		TS_ASSERT( m.has( 12 ) );
		TS_ASSERT( m[ 12 ] == 0 );
		TS_ASSERT( m.size() == 6 );

		// add 7th entry
		m.add( 22 );
		TS_ASSERT( m[ 22 ] == 0 );
		TS_ASSERT( m.size() == 7 );

		// test shrink
		m.shrink();
		TS_ASSERT( m.capacity() == m.size() );
	}


	/// @brief Swap test
	void test_SmallKeyMap_swap() {

		typedef  utility::keys::SmallKeyMap< int, int >  Map;
		Map m;
		m.assign( 3, 33 );
		m.assign( 6, 33 );
		m.assign( 9, 33 );
		Map w;
		w.assign( 3, 66 );
		w.assign( 6, 66 );
		w.assign( 9, 66 );
		w.assign( 12, 66 );
		w.assign( 15, 66 );
		w.assign( 18, 66 );

		// Swap maps
		swap( m, w );
		TS_ASSERT( m.size() == 6 );
		TS_ASSERT( w.size() == 3 );
		TS_ASSERT( m[ 3 ] == 66 );
		TS_ASSERT( w[ 9 ] == 33 );

		// Swap them back
		w.swap( m );
		TS_ASSERT( w.size() == 6 );
		TS_ASSERT( m.size() == 3 );
		TS_ASSERT( w[ 3 ] == 66 );
		TS_ASSERT( m[ 9 ] == 33 );
	}

};


