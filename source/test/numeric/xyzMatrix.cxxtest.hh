// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyzMatrix.cxxtest.hh
/// @brief  test suite for numeric::xyzMatrix
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.io.hh>
#include <iostream>


// --------------- Test Class --------------- //

class XYZMatrixTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.

	// three vectors used for initialization
	numeric::xyzVector_float v1;
	numeric::xyzVector_float v2;
	numeric::xyzVector_float v3;

	// two matrices that get filled with values 1-9
	numeric::xyzMatrix_float m1;  // filled by rows
	numeric::xyzMatrix_float m2;  // filled by columns

	// set tolerance for floating-point comparisons
	double delta_percent;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		using numeric::xyzVector_float;
		using numeric::xyzMatrix_float;

		// initialize vectors
		v1 = xyzVector_float( 1.0, 2.0, 3.0 );
		v2 = xyzVector_float( 4.0, 5.0, 6.0 );
		v3 = xyzVector_float( 7.0, 8.0, 9.0 );

		// initialize matrices
		m1 = xyzMatrix_float::rows( v1, v2, v3 );  // m1 has row-ordered values consecutive
		m2 = xyzMatrix_float::cols( v1, v2, v3 );  // m2 has column-ordered values consecutive

		delta_percent = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief Identity matrix constructor
	void test_xyzMatrix_ConstructIdentity() {

		using numeric::xyzMatrix_float;

		xyzMatrix_float m = xyzMatrix_float::identity();
		TS_ASSERT_EQUALS( m.xx(), 1.0 );
		TS_ASSERT_EQUALS( m.xy(), 0.0 );
		TS_ASSERT_EQUALS( m.xz(), 0.0 );

		TS_ASSERT_EQUALS( m.yx(), 0.0 );
		TS_ASSERT_EQUALS( m.yy(), 1.0 );
		TS_ASSERT_EQUALS( m.yz(), 0.0 );

		TS_ASSERT_EQUALS( m.zx(), 0.0 );
		TS_ASSERT_EQUALS( m.zy(), 0.0 );
		TS_ASSERT_EQUALS( m.zz(), 1.0 );
	}


	/// @brief Test construction from rows/cols
	void test_xyzMatrix_ConstructFromVectors() {

		using numeric::xyzMatrix_float;
		using numeric::xyzVector_float;

		// build a matrix from row vectors
		xyzMatrix_float m1( xyzMatrix_float::rows( v1, v2, v3 ) );
		TS_ASSERT_EQUALS( m1.xx(), 1.0 );
		TS_ASSERT_EQUALS( m1.xy(), 2.0 );
		TS_ASSERT_EQUALS( m1.xz(), 3.0 );
		TS_ASSERT_EQUALS( m1.yx(), 4.0 );
		TS_ASSERT_EQUALS( m1.yy(), 5.0 );
		TS_ASSERT_EQUALS( m1.yz(), 6.0 );
		TS_ASSERT_EQUALS( m1.zx(), 7.0 );
		TS_ASSERT_EQUALS( m1.zy(), 8.0 );
		TS_ASSERT_EQUALS( m1.zz(), 9.0 );

		// build a matrix from column vectors
		xyzMatrix_float m2( xyzMatrix_float::cols( v1, v2, v3 ) );
		TS_ASSERT_EQUALS( m2.xx(), 1.0 );
		TS_ASSERT_EQUALS( m2.xy(), 4.0 );
		TS_ASSERT_EQUALS( m2.xz(), 7.0 );
		TS_ASSERT_EQUALS( m2.yx(), 2.0 );
		TS_ASSERT_EQUALS( m2.yy(), 5.0 );
		TS_ASSERT_EQUALS( m2.yz(), 8.0 );
		TS_ASSERT_EQUALS( m2.zx(), 3.0 );
		TS_ASSERT_EQUALS( m2.zy(), 6.0 );
		TS_ASSERT_EQUALS( m2.zz(), 9.0 );

		// use the named row vector constructor
		xyzMatrix_float m3( xyzMatrix_float::rows_constructor( v1, v2, v3 ) );
		TS_ASSERT_EQUALS( m3, m1 );

		// use the named column vector constructor
		xyzMatrix_float m4( xyzMatrix_float::cols_constructor( v1, v2, v3 ) );
		TS_ASSERT_EQUALS( m4, m2 );
	}


	/// @brief Test construction from a list of values
	void test_xyzMatrix_ConstructFromValues() {

		using numeric::xyzMatrix_float;

		// put list of values in by rows
		xyzMatrix_float m1 = xyzMatrix_float::rows( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
		TS_ASSERT_EQUALS( m1, m1 );

		// put list of values in by columns
		xyzMatrix_float m2 = xyzMatrix_float::cols( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 );
		TS_ASSERT_EQUALS( m2, m2 );
	}


	/// @brief Test construction from a vector of contiguous values
	void test_xyzMatrix_ConstructFromContiguous() {

		using numeric::xyzMatrix_float;

		float v[ 10 ];
		for ( int i = 0; i <= 9; ++i ) {
			v[ i ] = i;
		}

		// test pointer to contiuguous values (row-ordered)
		xyzMatrix_float m1( xyzMatrix_float::rows( &v[ 1 ] ) );
		TS_ASSERT_EQUALS( m1, m1 );

		// test pointer to contiuguous values (row-ordered)
		xyzMatrix_float m2( xyzMatrix_float::cols( &v[ 1 ] ) );
		TS_ASSERT_EQUALS( m2, m2 );

		// test pointer to contiuguous row values
		xyzMatrix_float m3( xyzMatrix_float::rows( &v[ 1 ], &v[ 4 ], &v[ 7 ] ) );
		TS_ASSERT_EQUALS( m3, m1 );

		// test pointer to contiuguous column values
		xyzMatrix_float m4( xyzMatrix_float::cols( &v[ 1 ], &v[ 4 ], &v[ 7 ] ) );
		TS_ASSERT_EQUALS( m4, m2 );
	}


	/// @brief assignment using a single value
	void test_xyzMatrix_AssignFromSingleValue() {

		using numeric::xyzMatrix_float;

		xyzMatrix_float m1;
		m1 = 1.0f;
		TS_ASSERT_EQUALS( m1.xx(), 1.0 );
		TS_ASSERT_EQUALS( m1.xy(), 1.0 );
		TS_ASSERT_EQUALS( m1.xz(), 1.0 );
		TS_ASSERT_EQUALS( m1.yx(), 1.0 );
		TS_ASSERT_EQUALS( m1.yy(), 1.0 );
		TS_ASSERT_EQUALS( m1.yz(), 1.0 );
		TS_ASSERT_EQUALS( m1.zx(), 1.0 );
		TS_ASSERT_EQUALS( m1.zy(), 1.0 );
		TS_ASSERT_EQUALS( m1.zz(), 1.0 );
	}


	/// @brief assignment using contiguous values
	void test_xyzMatrix_AssignFromContiguous() {

		using numeric::xyzMatrix_float;

		float v[ 10 ];
		for ( int i = 0; i <= 9; ++i ) {
			v[ i ] = i;
		}

		// row-ordered assignment from contiguous values
		xyzMatrix_float m1;
		m1 = xyzMatrix_float::rows( &v[ 1 ] );
		TS_ASSERT_EQUALS( m1, m1 );

		// column-ordered assignment from contiguous values
		xyzMatrix_float m2;
		m2 = xyzMatrix_float::cols( &v[ 1 ] );
		TS_ASSERT_EQUALS( m2, m2 );
	}


	/// @brief assignment using rows/cols
	void test_xyzMatrix_AssignFromRowsCols() {

		using numeric::xyzMatrix_float;

		// Assign individual rows using row_*
		xyzMatrix_float m1;
		m1.row_x( v1 );
		m1.row_y( v2 );
		m1.row_z( v3 );
		TS_ASSERT_EQUALS( m1, m1 );

		// Assign individual columns using col_*
		xyzMatrix_float m2;
		m2.col_x( v1 );
		m2.col_y( v2 );
		m2.col_z( v3 );
		TS_ASSERT_EQUALS( m2, m2 );

		// assign individual rows using row index
		xyzMatrix_float m3;
		m3.row( 1, v1 );
		m3.row( 2, v2 );
		m3.row( 3, v3 );
		TS_ASSERT_EQUALS( m3, m1 );

		// assign individual columsn using column index
		xyzMatrix_float m4;
		m4.col( 1, v1 );
		m4.col( 2, v2 );
		m4.col( 3, v3 );
		TS_ASSERT_EQUALS( m4, m2 );
	}


	/// @brief assignment using rows/cols
	void test_xyzMatrix_Access() {

		using numeric::xyzMatrix_float;
		using numeric::xyzVector_float;

		// access rows using row_*
		xyzVector_float rx = m1.row_x();
		xyzVector_float ry = m1.row_y();
		xyzVector_float rz = m1.row_z();
		TS_ASSERT_EQUALS( rx, v1 );
		TS_ASSERT_EQUALS( ry, v2 );
		TS_ASSERT_EQUALS( rz, v3 );

		// access rows using row index
		xyzVector_float r1 = m1.row( 1 );
		xyzVector_float r2 = m1.row( 2 );
		xyzVector_float r3 = m1.row( 3 );
		TS_ASSERT_EQUALS( r1, v1 );
		TS_ASSERT_EQUALS( r2, v2 );
		TS_ASSERT_EQUALS( r3, v3 );

		// access columns using col_*
		xyzVector_float cx = m2.col_x();
		xyzVector_float cy = m2.col_y();
		xyzVector_float cz = m2.col_z();
		TS_ASSERT_EQUALS( cx, v1 );
		TS_ASSERT_EQUALS( cy, v2 );
		TS_ASSERT_EQUALS( cz, v3 );

		// access columns using column index
		xyzVector_float c1 = m2.col( 1 );
		xyzVector_float c2 = m2.col( 2 );
		xyzVector_float c3 = m2.col( 3 );
		TS_ASSERT_EQUALS( c1, v1 );
		TS_ASSERT_EQUALS( c2, v2 );
		TS_ASSERT_EQUALS( c3, v3 );
	}


	/// @brief test calculated properties
	void test_xyzMatrix_CalculatedProperties() {

		using numeric::xyzMatrix_float;

		xyzMatrix_float i = xyzMatrix_float::identity();

		// compute some determinants
		TS_ASSERT_DELTA( i.det(), 1.0f, delta_percent );
		TS_ASSERT_DELTA( m1.det(), 0.0f, delta_percent );

		// compute some traces
		TS_ASSERT_DELTA( i.trace(), 3.0f, delta_percent );
		TS_ASSERT_DELTA( m1.trace(), 15.0f, delta_percent );
	}


	/// @brief operations tests
	/// THIS TEST IS FAILING! LOOK INTO WHY!
	void test_xyzMatrix_Operations() {

		using numeric::xyzMatrix_float;
		using numeric::xyzVector_float;

		// multiply by a scalar
		xyzMatrix_float m1_doubled = xyzMatrix_float::rows( 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0 );
		TS_ASSERT_EQUALS( 2 * m1, m1_doubled );
		TS_ASSERT_EQUALS( m1 * 2, m1_doubled );

		// matrix/matrix product
		xyzMatrix_float m = m1;
		xyzMatrix_float m3 = xyzMatrix_float::rows( 1.0, 5.0, 2.0, 6.0, 8.0, 9.0, 4.0, 7.0, 3.0 );
		m.left_multiply_by( m3 );
		TS_ASSERT_DELTA( m.xx(),  35.0f, delta_percent );
		TS_ASSERT_DELTA( m.xy(),  43.0f, delta_percent );
		TS_ASSERT_DELTA( m.xz(),  51.0f, delta_percent );
		TS_ASSERT_DELTA( m.yx(), 101.0f, delta_percent );
		TS_ASSERT_DELTA( m.yy(), 124.0f, delta_percent );
		TS_ASSERT_DELTA( m.yz(), 147.0f, delta_percent );
		TS_ASSERT_DELTA( m.zx(),  53.0f, delta_percent );
		TS_ASSERT_DELTA( m.zy(),  67.0f, delta_percent );
		TS_ASSERT_DELTA( m.zz(),  81.0f, delta_percent );
	}

};


