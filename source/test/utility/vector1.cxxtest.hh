// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/vector1.cxxtest.hh
/// @brief  vector1.cxxtest: test suite for utility::vector1
/// @author Ron Jacak (ron.jacak@gmail.com)


// Testing headers
#include <cxxtest/TestSuite.h>
#include <utility/vector1.hh>

// Utility header
#include <utility/excn/EXCN_Base.hh>
#include <utility/backtrace.hh>

// C++ headers
#include <vector>

//namespace test {
//namespace utility {

// --- Types
typedef  std::vector< int >  std_vector_int;

class Vector1_Tests : public CxxTest::TestSuite {

	public:

		void test_vector1_constructor() {

			unsigned int vector_size = 5;
			utility::vector1_int v( vector_size, 99 );

			TS_ASSERT_EQUALS(v.size(), vector_size);
			TS_ASSERT_EQUALS(v[1], 99);
			TS_ASSERT_EQUALS(v[2], 99);
			TS_ASSERT_EQUALS(v[3], 99);
			TS_ASSERT_EQUALS(v[4], 99);
			TS_ASSERT_EQUALS(v[5], 99);

		}

		void test_vector1_assignment() {

			utility::vector1_int v( 3, 99 );
			utility::vector1_int w( 3, 88 );

			v = w;
			TS_ASSERT_EQUALS( v[1], 88 );
			TS_ASSERT_EQUALS( v[2], 88 );
			TS_ASSERT_EQUALS( v[3], 88 );
		}

		void test_vector1_copy() {

			utility::vector1_int v( 3, 99 );
			utility::vector1_int w( v );

			TS_ASSERT_EQUALS( w[1], 99 );
			TS_ASSERT_EQUALS( w[2], 99 );
			TS_ASSERT_EQUALS( w[3], 99 );

		}

		void test_vector1_compare_to_std() {

			utility::vector1_int v( 3, 99 );  // our version
			std_vector_int  w( 3, 99 );  // the std::vector version
			TS_ASSERT_EQUALS(v, w);
			TS_ASSERT( v == w );
			TS_ASSERT( !( v != w ) );
			TS_ASSERT( !( v < w ) );
			TS_ASSERT( v <= w );
			TS_ASSERT( !( v > w ) );
			TS_ASSERT( v >= w );
			TS_ASSERT( w == v );
			TS_ASSERT( !( w != v ) );
			TS_ASSERT( !( w < v ) );
			TS_ASSERT( w <= v );
			TS_ASSERT( !( w > v ) );
			TS_ASSERT( w >= v );

		}

		void test_vector1_swap() {

			// set up first vector
			utility::vector1_int v( 3 );
			v[1] = 1; v[2] = 2; v[3] = 3;
			TS_ASSERT_EQUALS(v[1], 1);
			TS_ASSERT_EQUALS(v[2], 2);
			TS_ASSERT_EQUALS(v[3], 3);

			// set up second vector and verify equal
			utility::vector1_int w( v );
			TS_ASSERT_EQUALS(v, w);

			// alter second vector and verify not equal
			w[1] = 5;
			TS_ASSERT_EQUALS(w[1], 5);
			TS_ASSERT_EQUALS(w[2], 2);
			TS_ASSERT_EQUALS(w[3], 3);
			TS_ASSERT( v != w );
			TS_ASSERT( w != v );

			// make fixed copies of each (for reference)
			utility::vector1_int const V( v );
			utility::vector1_int const W( w );
			TS_ASSERT_EQUALS(v, V);
			TS_ASSERT_EQUALS(w, W);

			// verify our utility::swap  <- depricated since gcc 4.2.*
			swap( v, w );
			TS_ASSERT_EQUALS(v, W);
			TS_ASSERT_EQUALS(w, V);

			// verify std::swap
			std::swap( v, w );
			TS_ASSERT_EQUALS(v, V);
			TS_ASSERT_EQUALS(w, W);

			// verify swap as method call
			v.swap( w );
			TS_ASSERT_EQUALS(v, W);
			TS_ASSERT_EQUALS(w, V);

			// Let C++ pick best swap match from std or utility
			// (This one might not actually test Koenig lookup, according to Ion?)
			using namespace std;
			swap( v, w );
			TS_ASSERT_EQUALS(v, V);
			TS_ASSERT_EQUALS(w, W);

		}


		/// @brief  Test assorted utility methods.
		/// @author Labonte <JWLabonte@jhu.edu>
		void test_vector1_utility_methods() {

			// Set up vectors.
			utility::vector1_int v( 3 );
			v[ 1 ] = 1; v[ 2 ] = 2; v[ 3 ] = 3;
			utility::vector1_int w( 3 );
			w[ 1 ] = 4; w[ 2 ] = 5; w[ 3 ] = 6;

			// Test append().
			v.append( w );
			TS_ASSERT_EQUALS( v.size(), 6 );
			TS_ASSERT_EQUALS( v[ 5 ], 5 );

			// Test contains().
			TS_ASSERT( v.contains( 4 ) );
			TS_ASSERT( ! v.contains( 0 ) );

			// Test index_of().
			TS_ASSERT_EQUALS( v.index_of( 1 ), 1 );
			TS_ASSERT_EQUALS( v.index_of( 6 ), 6 );
			try {
				set_throw_on_next_assertion_failure();
				v.index_of( 7 );  // This should force an exit.
				TS_ASSERT( false );  // Exception was not thrown!
			} catch ( utility::excn::EXCN_Base const & e) {
				std::string expected( "ERROR: vectorL:index_of: element not found\n" );
				TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ), expected.size() ), expected );
			}
		}


		/// @brief Test of bounds checking
		/// @note  This test is not used yet!  With the existing testing system,
		///        it would cause an assert() to fail and abort testing.  But perhaps
		///        it could be included in a future version that can catch runtime failures.
		//void test_vector1_bounds_check() {
		//	utility::vector1_int v( 5, 99 );
		//	TS_ASSERT_EQUALS(v[1000], 99); // should fail -- index out of bounds
		//}


};

//} // namespace utility
//} // namespace test

