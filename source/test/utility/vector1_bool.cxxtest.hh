// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/vector1_bool.cxxtest.hh
/// @brief  vector1_bool.cxxtest: test suite for utility::vector1_bool
/// @author Ron Jacak (ron.jacak@gmail.com)


// Testing headers
#include <cxxtest/TestSuite.h>
#include <utility/vector1_bool.hh>

// C++ headers
#include <vector>

// --- Types
typedef  std::vector< bool >  std_vector_bool;


class Vector1BoolTests : public CxxTest::TestSuite {

	public:

		void test_vector1_bool_constructor() {
			unsigned int vector_size = 5;
			utility::vector1_bool v( vector_size, false );

			TS_ASSERT_EQUALS(v.size(), vector_size);
			TS_ASSERT(v[1] == false);
			TS_ASSERT(v[2] == false);
			TS_ASSERT(v[3] == false);
			TS_ASSERT(v[4] == false);
			TS_ASSERT(v[5] == false);
		}


		/// @brief Assignment test
		/// VERIFY THIS ONE!
		void test_vector1_bool_assignment() {

			utility::vector1_bool v( 3, false );
			utility::vector1_bool w( 3, true );

			v = w;
			TS_ASSERT( v[1] == true );
			TS_ASSERT( v[2] == true );
			TS_ASSERT( v[3] == true );
		}


		/// @brief Copy constructor test
		/// VERIFY THIS ONE!
		void test_vector1_bool_copy() {

			utility::vector1_bool v( 3, false );
			utility::vector1_bool w( v );

			TS_ASSERT( w[1] == false );
			TS_ASSERT( w[2] == false );
			TS_ASSERT( w[3] == false );
		}


		/// @brief Test comparison to std::vector
		void test_vector1_bool_compare_to_std() {

			utility::vector1_bool v( 3, false );
			std_vector_bool  w( 3, false );

			TS_ASSERT(v == w);
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
		void test_vector1_bool_swap() {

			// set up first vector
			utility::vector1_bool v( 3 );  v[1] = true; v[2] = false; v[3] = true;
			TS_ASSERT(v[1] == true);
			TS_ASSERT(v[2] == false);
			TS_ASSERT(v[3] == true);

			// set up second vector and verify equal
			utility::vector1_bool w( v );
			TS_ASSERT_EQUALS(v, w);

			// alter second vector and verify not equal
			w[2] = true;
			TS_ASSERT(w[1] == true);
			TS_ASSERT(w[2] == true);
			TS_ASSERT(w[3] == true);
			TS_ASSERT( v != w );
			TS_ASSERT( w != v );

			// make fixed copies of each (for reference)
			utility::vector1_bool const V( v );
			utility::vector1_bool const W( w );
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


		/// @brief Flip test (specific to bool vectors)
		void test_vector1_bool_flip() {

			utility::vector1_bool v( 3 );  v[1] = true; v[2] = false; v[3] = true;

			// check starting condition
			TS_ASSERT(v[1] == true);
			TS_ASSERT(v[2] == false);
			TS_ASSERT(v[3] == true);

			// flip the whole vector
			v.flip();
			TS_ASSERT(v[1] == false);
			TS_ASSERT(v[2] == true);
			TS_ASSERT(v[3] == false);

			// flip a single element
			v[1].flip();
			TS_ASSERT(v[1] == true);
			TS_ASSERT(v[2] == true);
			TS_ASSERT(v[3] == false);

			// assign from an element
			v[3] = v[2];
			TS_ASSERT(v[1] == true);
			TS_ASSERT(v[2] == true);
			TS_ASSERT(v[3] == true);
		}


		/// @brief Test of bounds checking
		/// @note  This test is not used yet!  With the existing testing system,
		///        it would cause an assert() to fail and abort testing.  But perhaps
		///        it could be included in a future version that can catch runtime failures.
		//void test_vector1_bool_bounds_check() {
		//	utility::vector1_bool v( 5, false );
		//	TS_ASSERT_EQUALS(v[1000], false); // should fail -- index out of bounds
		//}


};

