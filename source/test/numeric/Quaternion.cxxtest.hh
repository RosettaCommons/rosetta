// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/Quaterion.cxxtest.hh
/// @brief  test suite for numeric::Quaternion
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/Quaternion.hh>
#include <numeric/Quaternion.io.hh>


typedef ::numeric::Quaternion_double QD;

// --------------- Test Class --------------- //

class QuaternionTests : public CxxTest::TestSuite {

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
	/// @brief Constructors
	void test_Quaternion_constructors() {

		QD i( QD::identity() );
		QD q( 1.0, 0.0, 0.0, 0.0 );

		// test various identity quaternions
		TS_ASSERT_EQUALS( q, i );
		TS_ASSERT_EQUALS( q, QD::I() );
		TS_ASSERT_EQUALS( q.magnitude(), 1.0 );
		TS_ASSERT_EQUALS( q.w(), 1.0 );
		TS_ASSERT_EQUALS( q.x(), 0.0 );
		TS_ASSERT_EQUALS( q.y(), 0.0 );
		TS_ASSERT_EQUALS( q.z(), 0.0 );
		TS_ASSERT_EQUALS( q.w_squared(), 1.0 );
		TS_ASSERT_EQUALS( q.x_squared(), 0.0 );
		TS_ASSERT_EQUALS( q.y_squared(), 0.0 );
		TS_ASSERT_EQUALS( q.z_squared(), 0.0 );

		// applying identity to itself should stay identity
		q.apply( q );
		TS_ASSERT_EQUALS( q, i );
		TS_ASSERT_EQUALS( q, QD::I() );
		TS_ASSERT_EQUALS( q.magnitude(), 1.0 );
	}

	/// @brief test normalization
	void test_Quaternion_normalize() {

		QD q( QD::identity() );
		const_cast< double & >( q.x() ) = .02; // Intentionally make non-unit
		q.normalize_if_needed( .1 ); // Shouldn't do anything
		TS_ASSERT( q.magnitude_squared() >= 1.000399 );
		q.normalize_if_needed( .0001 ); // Should normalize
		TS_ASSERT( q.magnitude_squared_error() < 1.0E-15 );
		TS_ASSERT_DELTA( q.magnitude(), 1.0, delta_percent );
	}

	/// @brief test conjugates
	void test_Quaternion_conjugate() {

		QD q1( 0.0, 1.0, 0.0, 0.0 );
		QD q2( q1.conjugated() );
		TS_ASSERT_EQUALS( q1 * q2, QD::I() );
		q2.conjugate();
		TS_ASSERT_EQUALS( q2, q1 );
	}

};


