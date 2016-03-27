// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vectorL.cxxtest.hh
/// @brief  vectorL.cxxtest: test suite for utility::vectorL
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <utility/vectorL.hh>
#include <cxxtest/TestSuite.h>

// C++ headers
#include <vector>

// --- Types
typedef  std::vector<int>  std_vector_int;


class VectorLTests : public CxxTest::TestSuite {

	public:

		/// @brief Size + Value Constructor Test (starting at -2)
		void test_vectorL_constructor_minus2() {
			unsigned int vector_size = 5;
			utility::vectorL< -2, int > v( vector_size, 99 );

			TS_ASSERT_EQUALS(v.size(), vector_size);
			TS_ASSERT_EQUALS(v[-2], 99);
			TS_ASSERT_EQUALS(v[-1], 99);
			TS_ASSERT_EQUALS(v[0], 99);
			TS_ASSERT_EQUALS(v[1], 99);
			TS_ASSERT_EQUALS(v[2], 99);
		}

		/// @brief Size + Value Constructor Test (starting at +2)
		void test_vectorL_constructor_plus2() {
			unsigned int vector_size = 5;
			utility::vectorL< 2, int > v( vector_size, 99 );

			TS_ASSERT_EQUALS(v.size(), vector_size);
			TS_ASSERT_EQUALS(v[2], 99);
			TS_ASSERT_EQUALS(v[3], 99);
			TS_ASSERT_EQUALS(v[4], 99);
			TS_ASSERT_EQUALS(v[5], 99);
			TS_ASSERT_EQUALS(v[6], 99);
		}


		/// @brief Test the has() method
		void test_vectorL_has_index() {
			unsigned int vector_size = 5;

			// test with bounds from -2 to 2
			utility::vectorL< -2, int >  v( vector_size, 99 );
			TS_ASSERT( v.has( -2 ) );
			TS_ASSERT( v.has( 2 ) );
			TS_ASSERT( ! v.has( -3 ) );
			TS_ASSERT( ! v.has( 3 ) );

			// test with bounds from +2 to +6
			utility::vectorL<  2, int >  w( vector_size, 99 );
			TS_ASSERT( w.has( 2 ) );
			TS_ASSERT( w.has( 6 ) );
			TS_ASSERT( ! w.has( 1 ) );
			TS_ASSERT( ! w.has( 7 ) );
		}


		/// @brief Assignment test
		/// VERIFY THIS ONE!
		void test_vectorL_assignment() {
			typedef  utility::vectorL< 0, int >  Vector;
			Vector v( 3, 99 );
			Vector w( 3, 88 );

			v = w;
			TS_ASSERT_EQUALS( v[0], 88 );
			TS_ASSERT_EQUALS( v[1], 88 );
			TS_ASSERT_EQUALS( v[2], 88 );
		}


		/// @brief Copy constructor test
		/// VERIFY THIS ONE!
		void test_vectorL_copy() {
			typedef  utility::vectorL< 0, int >  Vector;
			Vector v( 3, 99 );
			Vector w( v );

			TS_ASSERT_EQUALS( w[0], 99 );
			TS_ASSERT_EQUALS( w[1], 99 );
			TS_ASSERT_EQUALS( w[2], 99 );
		}


		/// @brief Test comparison to std::vector
		void test_vectorL_compare_to_std() {
			typedef  utility::vectorL< 0, int >  Vector;
			Vector v( 3, 99 );       // our version
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


		/// @brief Swap test
		void test_vectorL_swap() {
			// set up first vector
			typedef  utility::vectorL< 0, int >  Vector;
			Vector v( 3 );  v[0] = 1; v[1] = 2; v[2] = 3;
			TS_ASSERT_EQUALS(v[0], 1);
			TS_ASSERT_EQUALS(v[1], 2);
			TS_ASSERT_EQUALS(v[2], 3);

			// set up second vector and verify equal
			Vector w( v );
			TS_ASSERT_EQUALS(v, w);

			// alter second vector and verify not equal
			w[0] = 5;
			TS_ASSERT_EQUALS(w[0], 5);
			TS_ASSERT_EQUALS(w[1], 2);
			TS_ASSERT_EQUALS(w[2], 3);
			TS_ASSERT( v != w );
			TS_ASSERT( w != v );

			// make fixed copies of each (for reference)
			Vector const V( v );
			Vector const W( w );
			TS_ASSERT_EQUALS(v, V);
			TS_ASSERT_EQUALS(w, W);

			// verify our utility::swap  <- depricated since gcc 4.2.*
			utility::swap( v, w );
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
			//using namespace std;
			//swap( v, w );
			// <- this is no longer works due to swap call being ambiguous (after utility::swap was properly decalred as out-of-class friends), so instead use std::swap notation
			//TS_ASSERT_EQUALS(v, V);
			//TS_ASSERT_EQUALS(w, W);
		}


		/// @brief Shrink test
		void test_vectorL_shrink() {
			typedef  utility::vectorL< 0, int >  Vector;
			Vector v( 222, 33 );
			v.push_back( 44 ); // Generally this will cause capacity > size
			v.shrink();
			TS_ASSERT_EQUALS( v.size(), v.capacity() );
		}


		/// @brief Test of bounds checking
		/// @note  This test is not used yet!  With the existing testing system,
		///        it would cause an assert() to fail and abort testing.  But perhaps
		///        it could be included in a future version that can catch runtime failures.
		//void test_vectorL_bounds_check() {
		//	typedef utility::vectorL< 0, int > Vector;
		//	Vector v( 5, 99 );
		//	TS_ASSERT_EQUALS(v[1000], 99); // should fail -- index out of bounds
		//}

};  // class VectorLTests
