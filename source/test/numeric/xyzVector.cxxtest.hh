// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyzVector.cxxtest.hh
/// @brief  test suite for numeric::xyzVector
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyz.json.hh>

#include <boost/functional/hash.hpp>

#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/json_spirit/json_spirit_writer.h>

#include <sstream>
#include <iostream>

// --------------- Test Class --------------- //

class XYZVectorTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	double delta_percent;         // percentage difference for floating-point comparisons in TS_ASSERT_DELTA

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		delta_percent = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	/// @brief test construction from 1D array
	void test_xyzVector_ConstructFromSingleValue() {

		using numeric::xyzVector_float;

		xyzVector_float v( 15.0 );
		TS_ASSERT_EQUALS( v.x(), 15.0f );
		TS_ASSERT_EQUALS( v.y(), 15.0f );
		TS_ASSERT_EQUALS( v.z(), 15.0f );
	}


	/// @brief test construction from 1D array
	void test_xyzVector_ConstructFrom1DArrayPointer() {

		using numeric::xyzVector_float;

		// set up a 1D array
		float v[ 10 ];
		for ( int i = 0; i < 10; ++i ) {
			v[ i ] = i;
		}

		// initialize a vector from one spot
		xyzVector_float u( &v[ 1 ] );
		TS_ASSERT_EQUALS( u.x(), 1.0 );
		TS_ASSERT_EQUALS( u.y(), 2.0 );
		TS_ASSERT_EQUALS( u.z(), 3.0 );

		// initialize a vector from another spot
		xyzVector_float w( &v[ 4 ] );
		TS_ASSERT_EQUALS( w.x(), 4.0 );
		TS_ASSERT_EQUALS( w.y(), 5.0 );
		TS_ASSERT_EQUALS( w.z(), 6.0 );

		// Tests if constructor catches double to float conversion.
		// . caught at compilation time:
		// . no matching function for call to numeric::xyzVector<float>::xyzVector(double*)
		//FArray1D_double a( 3, 0.0 );
		//xyzVector_float b( &a( 1 ) );
	}


	/// @brief test construction from 2D array
	void test_xyzVector_ConstructFrom2DArrayPointer() {

		using numeric::xyzVector_float;
		// build a simple 2D array where the row values are sequential:
		//  0  1  2  3
		//  4  5  6  7
		//  8  9 10 11
		// 12 13 14 15
		float a [ 4 ] [ 4 ];
		for ( int i = 0; i < 4; ++i ) {
			for ( int j = 0; j < 4; ++j) {
				a[i][j] = 4*i + j;
			}
		}

		// initialize a vector from one spot
		xyzVector_float u( &a[0][0] );
		TS_ASSERT_EQUALS( u.x(), 0.0 );
		TS_ASSERT_EQUALS( u.y(), 1.0 );
		TS_ASSERT_EQUALS( u.z(), 2.0 );

		// initialize vector from another spot (across rows)
		xyzVector_float v( &a[2][3] );
		TS_ASSERT_EQUALS( v.x(), 11.0 );
		TS_ASSERT_EQUALS( v.y(), 12.0 );
		TS_ASSERT_EQUALS( v.z(), 13.0 );
	}


	/// @brief Length tests
	void test_xyzVector_Length() {

		using numeric::xyzVector_float;
		xyzVector_float v( 2.0, 1.0, -1.0 );

		// normalized vector should have length 1
		v.normalize();
		TS_ASSERT_DELTA( v.length(), 1.0f, delta_percent );
		// Note: First two parameters must be of same type for call to template function.

		// zeroed vector should have length 0
		v.zero();
		TS_ASSERT_DELTA( v.length(), 0.0f, delta_percent );
	}


	/// @brief Accessors tests
	void test_xyzVector_Accessors() {

		using numeric::xyzVector_float;
		float x = 1.0;
		float y = 2.0;
		float z = 3.0;

		xyzVector_float v( x, y, z );

		// test access to x/y/z components
		TS_ASSERT_EQUALS( v.x(), x );
		TS_ASSERT_EQUALS( v.y(), y );
		TS_ASSERT_EQUALS( v.z(), z );
	}


	/// @brief Relations tests
	void test_xyzVector_Relations() {

		using numeric::xyzVector_float;
		xyzVector_float v( 1.0, 2.0, 3.0 );
		xyzVector_float w( 1.0, 2.0, 3.0 );

		TS_ASSERT_EQUALS( v, w );

		// Reduce v and test inequality
		v -= 0.5;
		TS_ASSERT( v != w );
		TS_ASSERT( ! ( v == w ) );
		TS_ASSERT( v < w );
		TS_ASSERT( v <= w );

		// Increase v and test inequality
		v += 1.0;
		TS_ASSERT( v != w );
		TS_ASSERT( ! ( v == w ) );
		TS_ASSERT( v > w );
		TS_ASSERT( v >= w );

		// Test partial ordering: Set v.x to 0 but leave v.y > w.y and v.z > w.z so v and w are not orderable
		v.x( 0.0 );
		TS_ASSERT( v != w );
		TS_ASSERT( ! ( v == w ) );
		TS_ASSERT( ! ( v < w ) );
		TS_ASSERT( ! ( v <= w ) );
		TS_ASSERT( ! ( v > w ) );
		TS_ASSERT( ! ( v >= w ) );

		// Test length relations
		TS_ASSERT( ! equal_length( v, w ) );
		TS_ASSERT( ! v.equal_length(w) );
		TS_ASSERT( not_equal_length( v, w ) );
		TS_ASSERT( v.not_equal_length(w) );
		TS_ASSERT( v.longer(w) );
		TS_ASSERT( ! v.shorter(w) );
		TS_ASSERT( v.longer_or_equal(w) );
		TS_ASSERT( ! v.shorter_or_equal(w) );
	}


	/// @brief binary operator tests
	void test_xyzVector_BinaryOperations() {

		using numeric::xyzVector_float;

		xyzVector_float v( 1.0, 2.0, 3.0 );
		xyzVector_float w( 1.0, 2.0, 3.0 );

		// Check dot product of two equal vectors
		float x = dot( v, w );
		TS_ASSERT_DELTA( x, v.length_squared(), delta_percent ); // v == w here

		// Tweak the vectors and compute cross product
		v += 1.0; w -= 1.0;
		xyzVector_float t = cross( v, w );
		TS_ASSERT_DELTA( dot( t, v ), 0.0f, delta_percent ); // t and v are orthogonal
		TS_ASSERT_DELTA( dot( t, w ), 0.0f, delta_percent ); // t and w are orthogonal

		// Check midpoint (should match original vector)
		xyzVector_float original( 1.0, 2.0, 3.0 );
		xyzVector_float mid;
		mid = midpoint ( v, w );
		TS_ASSERT_DELTA( original.x(), mid.x(), delta_percent );
		TS_ASSERT_DELTA( original.y(), mid.y(), delta_percent );
		TS_ASSERT_DELTA( original.z(), mid.z(), delta_percent );
	}


	/// @brief boost::hash tests
	void test_xyzVector_HashOperations() {
		numeric::xyzVector_float v( 1.0, 2.0, 3.0 );
		numeric::xyzVector_float w( 1.0, 2.0, 3.0 );
		numeric::xyzVector_float x( 2.0, 2.0, 3.0 );


		boost::hash<numeric::xyzVector_float> xyz_hasher;

		TS_ASSERT(xyz_hasher(v) == xyz_hasher(w));
		TS_ASSERT(xyz_hasher(v) != xyz_hasher(x));

	}

	/// @brief serialization tests
	void test_xyzVector_Serialization() {
		numeric::xyzVector_float v( 1.0, 2.0, 3.0 );

		std::ostringstream ss;
		ss << utility::json_spirit::write(numeric::serialize(v));
		utility::json_spirit::mValue v_data;
		utility::json_spirit::read(ss.str(),v_data);
		numeric::xyzVector_float w;
		utility::json_spirit::mArray point_array(v_data.get_array());
		w = numeric::deserialize<float>(point_array);
		TS_ASSERT_DELTA(v.x(), w.x(),delta_percent);
		TS_ASSERT_DELTA(v.y(), w.y(),delta_percent);
		TS_ASSERT_DELTA(v.z(), w.z(),delta_percent);

	}

};

