// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/small_vector1.cxxtest.hh
/// @brief  boost::container::small_vector with 1-based indexing
/// @author Jack Maguire, jackmaguire1444@gmail.com


// Testing headers
#include <cxxtest/TestSuite.h>
#include <utility/small_vector1.hh>

// Utility header
#include <utility/excn/Exceptions.hh>
#include <utility/backtrace.hh>

// C++ headers
#include <vector>

//namespace test {
//namespace utility {

// --- Types
typedef  std::vector< int >  std_vector_int;

class SmallVector1_Tests : public CxxTest::TestSuite {

	public:

		void test_vector1_constructor() {

			unsigned int vector_size = 5;
			utility::small_vector1< int, 2 > v( vector_size, 99 );

			TS_ASSERT_EQUALS(v.size(), vector_size);
			TS_ASSERT_EQUALS(v[1], 99);
			TS_ASSERT_EQUALS(v[2], 99);
			TS_ASSERT_EQUALS(v[3], 99);
			TS_ASSERT_EQUALS(v[4], 99);
			TS_ASSERT_EQUALS(v[5], 99);

		}

		void test_vector1_assignment() {

			utility::small_vector1< int, 2 > v( 3, 99 );
			utility::small_vector1< int, 6 > w( 3, 88 );

			v = w;
			TS_ASSERT_EQUALS( v[1], 88 );
			TS_ASSERT_EQUALS( v[2], 88 );
			TS_ASSERT_EQUALS( v[3], 88 );
		}

		void test_vector1_copy() {

			utility::small_vector1< int, 2 > v( 3, 99 );
			utility::small_vector1< int, 2 > w( v );

			TS_ASSERT_EQUALS( w[1], 99 );
			TS_ASSERT_EQUALS( w[2], 99 );
			TS_ASSERT_EQUALS( w[3], 99 );

		}

		void test_vector1_swap() {

			// set up first vector
			utility::small_vector1< int, 2 > v( 3 );
			v[1] = 1; v[2] = 2; v[3] = 3;
			TS_ASSERT_EQUALS(v[1], 1);
			TS_ASSERT_EQUALS(v[2], 2);
			TS_ASSERT_EQUALS(v[3], 3);

			// set up second vector and verify equal
			utility::small_vector1< int, 2 > w( v );
			TS_ASSERT_EQUALS(v, w);

			// alter second vector and verify not equal
			w[1] = 5;
			TS_ASSERT_EQUALS(w[1], 5);
			TS_ASSERT_EQUALS(w[2], 2);
			TS_ASSERT_EQUALS(w[3], 3);
			TS_ASSERT( v != w );
			TS_ASSERT( w != v );

			// make fixed copies of each (for reference)
			utility::small_vector1< int, 2 > const V( v );
			utility::small_vector1< int, 2 > const W( w );
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

};

//} // namespace utility
//} // namespace test

