// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/ClassKeyMap.unit.cc
/// @brief  ClassKeyMap unit test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Package headers
#include <utility/keys/ClassKeyMap.test.hh>
#include <utility/keys/ClassKeyMap.hh>


namespace test {
namespace utility {
namespace keys {


// --- set up the test suite
TEST_SUITE_BEGIN(ClassKeyMap_Tests)
	TEST_SUITE_USES_CASE(test_ClassKeyMap_general)
	TEST_SUITE_USES_CASE(test_ClassKeyMap_swap)
TEST_SUITE_END


// --- set up the individual test cases


/// @brief General tests for int
TEST_CASE_BEGIN(test_ClassKeyMap_general)
{
	typedef  ::utility::keys::ClassKeyMap< int, int, float >  Map;

	// Create an instance and check starting conditions
	Map m;
	m.assign( 1, 123 ).assign( 2, 123 ).assign( 3, 123 );
	TEST_CHECK( m.size() == 3 );

	// Assign another key/value pair
	m.assign( 7, 77 );
	TEST_CHECK( m[ 7 ] == 77 );

	// Add a value with operator()
	m( 8 ) = 88;
	TEST_CHECK( m[ 8 ] == 88 );

	// Add a key
	m.add( 12 );
	TEST_CHECK( m.has( 12 ) );
	TEST_CHECK( m[ 12 ] == 0 );
	TEST_CHECK( m.size() == 6 );

	// Add another key
	m.add( 22 );
	TEST_CHECK( m.has( 22 ) );
	TEST_CHECK( m[ 22 ] == 0 );
	TEST_CHECK( m.size() == 7 );

	// Shrink to right-sized
	m.shrink();
	TEST_CHECK( m.capacity() == m.size() );
}
TEST_CASE_END


/// @brief Swap test
TEST_CASE_BEGIN(test_ClassKeyMap_swap)
{
	typedef  ::utility::keys::ClassKeyMap< int, int, double >  Map;

	// Create an instance and check starting conditions
	Map m;
	m.assign( 3, 33 ).assign( 6, 33 ).assign( 9, 33 );
	TEST_CHECK( m.size() == 3 );

	// Creat a second ClassKeyMap
	Map w;
	w.assign( 3, 66 ).assign( 6, 66 ).assign( 9, 66 ).assign( 12, 66 ).assign( 15, 66 ).assign( 18, 66 );
	TEST_CHECK( w.size() == 6 );

	// Swap the two ClassKeyMaps
	swap( m, w );
	TEST_CHECK( m.size() == 6 );
	TEST_CHECK( w.size() == 3 );
	TEST_CHECK( m[ 3 ] == 66 );
	TEST_CHECK( w[ 9 ] == 33 );

	// Swap the two ClassKeyMaps back
	w.swap( m );
	TEST_CHECK( w.size() == 6 );
	TEST_CHECK( m.size() == 3 );
	TEST_CHECK( w[ 3 ] == 66 );
	TEST_CHECK( m[ 9 ] == 33 );
}
TEST_CASE_END


} // namespace keys
} // namespace utility
} // namespace test
