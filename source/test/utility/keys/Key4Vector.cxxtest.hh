// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key4Vector.cxxtest.hh
/// @brief  Key4Vector unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/keys/Key4Vector.hh>
#include <utility/keys/Key4Tuple.hh>


class Key4VectorTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests
	void test_Key4Vector_general() {

		using namespace utility::keys;

		typedef  Key4Vector< int >  Key4Vector_int;

		{ // Constructor and accessors
			Key4Vector_int key( 1, 2, 3, 4 );
			TS_ASSERT( key[ 0 ] == 1 );
			TS_ASSERT( key[ 1 ] == 2 );
			TS_ASSERT( key[ 2 ] == 3 );
			TS_ASSERT( key[ 3 ] == 4 );
			TS_ASSERT( key( 1 ) == 1 );
			TS_ASSERT( key( 2 ) == 2 );
			TS_ASSERT( key( 3 ) == 3 );
			TS_ASSERT( key( 4 ) == 4 );
			TS_ASSERT( key.key1() == 1 );
			TS_ASSERT( key.key2() == 2 );
			TS_ASSERT( key.key3() == 3 );
			TS_ASSERT( key.key4() == 4 );
		}

		{ // Copy constructor and comparisons
			Key4Vector_int keyA( 1, 2, 3, 4 ), keyB( keyA );
			TS_ASSERT( keyA == keyB );
			TS_ASSERT( ! ( keyA < keyB ) );
			TS_ASSERT( ! ( keyB < keyA ) );
		}

		{ // Default constructor, assignment, and comparisons
			Key4Vector_int keyA( 1, 2, 3, 4 ), keyB;
			keyB = keyA;
			TS_ASSERT( keyA == keyB );
			TS_ASSERT( ! ( keyA != keyB ) );
			TS_ASSERT( ! ( keyA < keyB ) );
			TS_ASSERT( ! ( keyB < keyA ) );
		}

		{ // Default constructor, assignment, and comparisons
			Key4Vector_int keyA( 1, 2, 3, 4 ), keyB;
			keyB = keyA;
			TS_ASSERT( keyA == keyB );
			TS_ASSERT( ! ( keyA != keyB ) );
			TS_ASSERT( ! ( keyA < keyB ) );
			TS_ASSERT( ! ( keyB < keyA ) );
		}

		{ // Lexicographic ordering
			Key4Vector_int keyA( 1, 2, 3, 4 );
			Key4Vector_int keyB( 1, 2, 4, 4 );
			TS_ASSERT( ! ( keyA == keyB ) );
			TS_ASSERT( keyA != keyB );
			TS_ASSERT( keyA < keyB );
			TS_ASSERT( ! ( keyB < keyA ) );
		}

		{ // Lexicographic ordering
			Key4Vector_int keyA( 1, 2, 3, 4 );
			Key4Vector_int keyB( 2, 1, 2, 3 );
			TS_ASSERT( ! ( keyA == keyB ) );
			TS_ASSERT( keyA != keyB );
			TS_ASSERT( keyA < keyB );
			TS_ASSERT( ! ( keyB < keyA ) );
		}

	}

};


