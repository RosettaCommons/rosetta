// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/ClassKeyVector.cxxtest.hh
/// @brief  ClassKeyVector unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/keys/ClassKeyVector.hh>


class KeyVectorTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests for int
	void test_ClassKeyVector_int() {

		typedef  utility::keys::ClassKeyVector< int, int, int >  Vector;

		// create an instance and check starting conditions
		Vector v( 3, 123 );
		TS_ASSERT( v.size() == 3 );
		TS_ASSERT( ! v.has( 7 ) );

		// add a key/value assignment to the vector
		v.assign( 7, 77 );
		TS_ASSERT( v[ 7 ] == 77 );
		TS_ASSERT( v.has( 7 ) );

		// add a second key/value entry
		v( 8 ) = 88;
		TS_ASSERT( v[ 8 ] == 88 );
		TS_ASSERT( v.has( 8 ) );

		// activate a third entry (but keep default value)
		Vector::activate( 12 );
		TS_ASSERT( v[ 12 ] == 123 );

		// verify that the vector size is still only 3
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


	/// @brief General tests for long
	void test_ClassKeyVector_long() {

		typedef  utility::keys::ClassKeyVector< long, int, int >  Vector;

		// create an instance and check starting conditions
		Vector v( 3, 123 );
		TS_ASSERT( v.size() == 3 );

		// add a key/value assignment to the vector
		v.assign( 7l, 77 );
		TS_ASSERT( v[ 7l ] == 77 );

		// add a second key/value entry
		v( 8l ) = 88;

		// activate a third entry (but keep default value)
		Vector::activate( 12l );
		TS_ASSERT( v[ 12l ] == 123 );
		TS_ASSERT( v.size() == 3 );

		// add a fourth entry, which should expand the vector size
		v.add( 22l );
		TS_ASSERT( v.size() == 4 );

		// we expanded past original 3 entries, so 4th value should be zero
		TS_ASSERT( v[ 22l ] == 0 );

		// shrink the vector
		v.shrink();
		TS_ASSERT( v.capacity() == v.size() );
	}


	/// @brief Swap test
	void test_ClassKeyVector_swap() {

		typedef  utility::keys::ClassKeyVector< int, int, double >  Vector;
		Vector v( 3, 33 ), w( 6, 66 );
		Vector::activate( 3 );
		Vector::activate( 6 );
		Vector::activate( 9 );

		// Swap vectors
		swap( v, w );
		TS_ASSERT( v.size() == 6 );
		TS_ASSERT( w.size() == 3 );
		TS_ASSERT( v[ 3 ] == 66 );
		TS_ASSERT( w[ 9 ] == 33 );

		// Swap them back
		w.swap( v );
		TS_ASSERT( w.size() == 6 );
		TS_ASSERT( v.size() == 3 );
		TS_ASSERT( w[ 3 ] == 66 );
		TS_ASSERT( v[ 9 ] == 33 );
	}

};


