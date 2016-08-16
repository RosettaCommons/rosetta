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

// Utility headers
#include <utility/vector1.hh>

// Package headers
#include <numeric/PCA.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.io.hh>
#include <iostream>


// --------------- Test Class --------------- //

class PCATests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.

	// points to run PCA on
	utility::vector1<numeric::xyzVector_float> points;
	
	double delta_percent;         // percentage difference for floating-point comparisons in TS_ASSERT_DELTA

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		using numeric::xyzVector_float;

		// initialize points
		points.push_back(xyzVector_float(1.0,2.0,3.0));
		points.push_back(xyzVector_float(3.0,2.0,1.0));
		points.push_back(xyzVector_float(4.0,5.0,6.0));

		delta_percent = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief Identity matrix constructor
	void test_pca() {

		using numeric::xyzMatrix_float;

		xyzMatrix_float pca_results =
			principal_components(points);

		TS_ASSERT_DELTA( pca_results.xx(), 0.31775412, delta_percent );
		TS_ASSERT_DELTA( pca_results.yx(), 0.54635078, delta_percent );
		TS_ASSERT_DELTA( pca_results.zx(), 0.77494085, delta_percent );

		TS_ASSERT_DELTA( pca_results.xy(), 0.85530704, delta_percent );
		TS_ASSERT_DELTA( pca_results.yy(), 0.18759693, delta_percent );
		TS_ASSERT_DELTA( pca_results.zy(), -0.48296729, delta_percent );

		TS_ASSERT_DELTA( pca_results.xz(), -0.40924600, delta_percent );
		TS_ASSERT_DELTA( pca_results.yz(), 0.81627715, delta_percent );
		TS_ASSERT_DELTA( pca_results.zz(), -0.40768793, delta_percent );

	}

	void test_first_pc() {

		using numeric::xyzVector_float;

		xyzVector_float first_pc =
			first_principal_component(points);

		TS_ASSERT_DELTA( first_pc.x(), 0.31775412, delta_percent );
		TS_ASSERT_DELTA( first_pc.y(), 0.54635078, delta_percent );
		TS_ASSERT_DELTA( first_pc.z(), 0.77494085, delta_percent );

	}

};


