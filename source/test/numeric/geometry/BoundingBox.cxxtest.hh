// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/BoundingBox.cxxtest.hh
/// @brief  test suite for numeric::geometry::BoundingBox
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/xyzVector.hh>

// --------------- Test Class --------------- //

class BoundingBoxTests : public CxxTest::TestSuite {

public:

	// typedefs
	typedef numeric::xyzVector< double > Point;
	typedef numeric::geometry::BoundingBox< Point > BoundingBox;

	// shared data
	BoundingBox bb;
	Point bb1_min;
	Point bb1_max;

	// shared initialization
	void setUp() {
		// setup dummy bounding box
		bb = BoundingBox( Point( 0.0, 0.0, 0.0 ) );
		bb1_min = Point( -5.67, -7.7, -4.33 );
		bb1_max = Point( 8.9, 3.12, 5.55 );

		bb.add( bb1_min );
		bb.add( bb1_max );
	}

	// shared finalization
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	/// @brief test lower stays the same upon addition of point
	void test_BoundingBox_LowerRemains() {
		Point p( 0.34, -7.0, -4.33 );
		bb.add( p );
		TS_ASSERT_EQUALS( bb.lower().x(), bb1_min.x() );
		TS_ASSERT_EQUALS( bb.lower().y(), bb1_min.y() );
		TS_ASSERT_EQUALS( bb.lower().z(), bb1_min.z() );
	}


	/// @brief test upper stays the same upon addition of point
	void test_BoundingBox_UpperRemains() {
		Point p( 8.9, -7.0, 1.2 );
		bb.add( p );
		TS_ASSERT_EQUALS( bb.upper().x(), bb1_max.x() );
		TS_ASSERT_EQUALS( bb.upper().y(), bb1_max.y() );
		TS_ASSERT_EQUALS( bb.upper().z(), bb1_max.z() );
	}


	/// @brief test point containment
	void test_BoundingBox_ContainsPoint() {
		Point p( 1.1, -1.2, 0.55 );
		TS_ASSERT( bb.contains( p ) );
	}


	/// @brief test point non-containment
	void test_BoundingBox_DoesNotContainPoint() {
		Point p( 10.0, -10.0, 10.0 );
		TS_ASSERT( !bb.contains( p ) );
	}


	/// @brief test box intersection (partial containment)
	void test_BoundingBox_IntersectsPartialContainment() {
		BoundingBox bb2( Point( 0.0, 0.0, 0.0 ) );
		bb2.add( Point( -10.0, -12.2, -4.0 ) );
		TS_ASSERT( bb.intersects( bb2 ) );
	}


	/// @brief test box intersection (full containment)
	void tests_BoundingBox_IntersectsFullContainment() {
		BoundingBox bb2( Point( 0.0, 0.0, 0.0 ) );
		bb2.add( Point( -1.0, -1.0, -1.0 ) );
		bb2.add( Point( 1.0, 1.0, 1.0 ) );
		TS_ASSERT( bb.intersects( bb2 ) );
	}


	/// @brief test box non-intersection
	void test_BoundingBox_DoesNotIntersect() {
		BoundingBox bb2( Point( 10.0, 10.0, 10.0 ) );
		bb2.add( Point( 15.0, 15.0, 15.0 ) );
		TS_ASSERT( !bb.intersects( bb2 ) );
	}

};
