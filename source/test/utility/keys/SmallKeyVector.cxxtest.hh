// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/SmallKeyVector.cxxtest.hh
/// @brief  SmallKeyVector unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/keys/SmallKeyVector.hh>


class SmallKeyVectorTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests
	void test_SmallKeyVector_general() {

		typedef  utility::keys::SmallKeyVector< int, int >  Vector;

		// create an instance and check starting conditions
		Vector v( 3, 123 );
		TS_ASSERT( v.size() == 3 );

		// add a key/value assignment to the vector
		v.assign( 7, 77 );
		TS_ASSERT( v[ 7 ] == 77 );

		// add a second key/value entry
		v( 8 ) = 88;

		// activate a third entry (but keep default value)
		v.add( 12 );
		TS_ASSERT( v[ 12 ] == 123 );
		TS_ASSERT( v.size() == 3 );

		// add a fourth entry, which should expand the vector size
		v.add( 22 );
		TS_ASSERT( v.size() == 4 );

		// we expanded past original 3 entries, so 4th value should be zero
		TS_ASSERT( v[ 22 ] == 0 );

		// shrink the vector
		v.shrink();
		TS_ASSERT( v.capacity() == v.size() );
	}


	/// @brief Swap test
	void test_SmallKeyVector_swap() {

		typedef  utility::keys::SmallKeyVector< int, int >  Vector;
		Vector v( 3, 33 ), w( 6, 66 );
		v.add( 3 );
		v.add( 6 );
		v.add( 9 );

		// Swap vectors
		swap( v, w );
		TS_ASSERT( v.size() == 6 );
		TS_ASSERT( w.size() == 3 );
		TS_ASSERT( v[ Vector::Index( 3u ) ] == 66 ); // 3u gets Index, not Key, accessor
		TS_ASSERT( w[ 9 ] == 33 );

		// Swap them back
		w.swap( v );
		TS_ASSERT( w.size() == 6 );
		TS_ASSERT( v.size() == 3 );
		TS_ASSERT( w[ Vector::Index( 3u ) ] == 66 ); // 3u gets Index, not Key, accessor
		TS_ASSERT( v[ 9 ] == 33 );
	}

};


