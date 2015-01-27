// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vector0.cxxtest.hh
/// @brief  vector0.test: test suite for utility::vector0
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/vector0.hh>

// C++ headers
#include <vector>


// --- Types
typedef  std::vector<int>  std_vector_int;


class Vector0Tests : public CxxTest::TestSuite {

	public:

		/// @brief Size + Value Constructor Test
		void test_vector0_constructor() {

			unsigned int vector_size = 5;
			utility::vector0_int v( vector_size, 99 );

			TS_ASSERT_EQUALS(v.size(), vector_size);
			TS_ASSERT_EQUALS(v[0], 99);
			TS_ASSERT_EQUALS(v[1], 99);
			TS_ASSERT_EQUALS(v[2], 99);
			TS_ASSERT_EQUALS(v[3], 99);
			TS_ASSERT_EQUALS(v[4], 99);
		}


		/// @brief Assignment test
		void test_vector0_assignment() {

			utility::vector0_int v( 3, 99 );
			utility::vector0_int w( 3, 88 );

			v = w;
			TS_ASSERT_EQUALS( v[0], 88 );
			TS_ASSERT_EQUALS( v[1], 88 );
			TS_ASSERT_EQUALS( v[2], 88 );
		}


		/// @brief Copy constructor test
		void test_vector0_copy() {

			utility::vector0_int v( 3, 99 );
			utility::vector0_int w( v );

			TS_ASSERT_EQUALS( w[0], 99 );
			TS_ASSERT_EQUALS( w[1], 99 );
			TS_ASSERT_EQUALS( w[2], 99 );
		}


		/// @brief Test comparison to std::vector
		void test_vector0_compare_to_std() {

			utility::vector0_int v( 3, 99 );
			std_vector_int  w( 3, 99 );

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


		/// @brief Swap test
		void test_vector0_swap() {

			// set up first vector
			utility::vector0_int v( 3 );  v[0] = 1; v[1] = 2; v[2] = 3;
			TS_ASSERT_EQUALS(v[0], 1);
			TS_ASSERT_EQUALS(v[1], 2);
			TS_ASSERT_EQUALS(v[2], 3);

			// set up second vector and verify equal
			utility::vector0_int w( v );
			TS_ASSERT_EQUALS(v, w);

			// alter second vector and verify not equal
			w[0] = 5;
			TS_ASSERT_EQUALS(w[0], 5);
			TS_ASSERT_EQUALS(w[1], 2);
			TS_ASSERT_EQUALS(w[2], 3);
			TS_ASSERT( v != w );
			TS_ASSERT( w != v );

			// make fixed copies of each (for reference)
			utility::vector0_int const V( v );
			utility::vector0_int const W( w );
			TS_ASSERT_EQUALS(v, V);
			TS_ASSERT_EQUALS(w, W);

			// verify our swap
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
		void test_vector0_utility_methods() {

			// Set up vectors.
			utility::vector0_int v( 3 );
			v[ 0 ] = 1; v[ 1 ] = 2; v[ 2 ] = 3;
			utility::vector0_int w( 3 );
			w[ 0 ] = 4; w[ 1 ] = 5; w[ 2 ] = 6;
			
			// Test append().
			v.append( w );
			TS_ASSERT_EQUALS( v.size(), 6 );
			TS_ASSERT_EQUALS( v[ 4 ], 5 );
			
			// Test contains().
			TS_ASSERT( v.contains( 4 ) );
			TS_ASSERT( ! v.contains( 0 ) );
			
			// Test index_of().
			TS_ASSERT_EQUALS( v.index_of( 1 ), 0 );
			TS_ASSERT_EQUALS( v.index_of( 6 ), 5 );
			TS_TRACE( "An out-of-bounds error message should follow:" );
			try {
				v.index_of( 7 );  // This should force an exit.
				TS_ASSERT( false );  // Exception was not thrown!
			} catch ( utility::excn::EXCN_Base const & e) {
				TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ) ),
						"ERROR: vectorL:index_of: element not found\n\n" );
				TS_TRACE( "The above error message was expected." );
			}
		}

		/// @brief Test of bounds checking
		/// @note  This test is not used yet!  With the existing testing system,
		///        it would cause an assert() to fail and abort testing.  But perhaps
		///        it could be included in a future version that can catch runtime failures.
		//void test_vector0_bounds_check() {
		//	::utility::vector0_int v( 5, 99 );
		//	TS_ASSERT_EQUALS(v[1000], 99); // should fail -- index out of bounds
		//}


}; // class Vector0_Tests

