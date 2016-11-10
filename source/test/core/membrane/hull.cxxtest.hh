// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/membrane/hull.cxxtest.hh
/// @brief   Unit test for hull util functions
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/membrane/hull.hh>

// Package Headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

static THREAD_LOCAL basic::Tracer TR("core.membrane.hull.cxxtest");

using namespace core::membrane;

class HullTest : public CxxTest::TestSuite {

public: // test functions

	/// Test Setup Functions ////////

	/// @brief Setup Test
	void setUp(){

		// Initialize
		core_init();
	}

	/// @brief Standard Tear Down
	void tearDown(){}

	///// Test Methods /////////////

	////////////////////////////////////////////////////////////////////////////////

	// does point q to infinity intersect with p1p2 in 2D?
	void test_intersect() {

		TR << "test intersect" << std::endl;

		core::Vector p1( 1, 1, 1 );
		core::Vector p2( 7, 7, 7 );
		core::Vector q1( 2, 5, 3 );
		core::Vector q2( 9, 9, 9 );

		TS_ASSERT_EQUALS( intersect( p1, p2, q1 ), true );
		TS_ASSERT_EQUALS( intersect( p1, p2, q2 ), false );

		core::Vector p3( 7, 0, 10);
		core::Vector p4( 7, 7, 7 );
		core::Vector q3( -4, 3, 8 );

		TS_ASSERT_EQUALS( intersect( p3, p4, q3 ), true );

	}
	
	// does q lie on segment p1p2 in 2D?
	void test_on_segment() {
		
		TR << "test on segment" << std::endl;
		
		core::Vector p1( 1, 1, 1 );
		core::Vector p2( 7, 7, 7 );
		core::Vector q1( 2, 2, 2 );
		core::Vector q2( 3, 4, 5 );

		TS_ASSERT_EQUALS( on_segment( p1, p2, q1 ), true );
		TS_ASSERT_EQUALS( on_segment( p1, p2, q2 ), false );

		core::Vector p3( 7, 0, 10);
		core::Vector p4( 7, 7, 7 );
		core::Vector q4( 7, 5, 12 );
		TS_ASSERT_EQUALS( on_segment( p3, p4, q4 ), true );

	}

	// does point q lie inside 2D polygon defined by points in map
	void test_inside_polygon() {

		TR << "test inside polygon" << std::endl;
		
		core::Vector p1( 1, 1, 1 );
		core::Vector p2( 7, 0, 10);
		core::Vector p3( 7, 7, 7 );
		core::Vector p4( 2, 8, 2 );
		
		std::map< core::Size, core::Vector > coords;
		coords[ 1 ] = p1;
		coords[ 2 ] = p2;
		coords[ 3 ] = p3;
		coords[ 4 ] = p4;
		
		utility::vector1< core::Size > polygon;
		polygon.push_back( 2 );
		polygon.push_back( 3 );
		polygon.push_back( 4 );
		polygon.push_back( 1 );
		
		core::Vector q1( 3, 4, 5 );
		core::Vector q2( 5, 20, 4 );
		core::Vector q3( -4, 3, 8 );
		core::Vector q4( 7, 5, 12 );
		
		TS_ASSERT_EQUALS( inside_polygon( coords, polygon, q1 ), "i" );
		TS_ASSERT_EQUALS( inside_polygon( coords, polygon, q2 ), "o" );
		TS_ASSERT_EQUALS( inside_polygon( coords, polygon, q3 ), "o" );
		TS_ASSERT_EQUALS( inside_polygon( coords, polygon, q4 ), "b" );

	}

	// test whether point q is to left of p1p2
	void test_to_left() {
		
		TR << "test to left" << std::endl;
		
		core::Vector p1( 1, 1, 1 );
		core::Vector p2( 7, 7, 10);
		
		core::Vector q1( 3, 4, 5 );
		core::Vector q2( 8, 3, 1 );

		TS_ASSERT_EQUALS( to_left( p1, p2, q1 ), true );
		TS_ASSERT_EQUALS( to_left( p1, p2, q2 ), false );
		TS_ASSERT_EQUALS( to_left( p2, p1, q1 ), true );
		TS_ASSERT_EQUALS( to_left( p2, p1, q2 ), false );

	}
	
	// test whether q is clockwise from p1p2
	void test_clockwise() {
		
		TR << "test clockwise" << std::endl;
		
		core::Vector p1( 1, 1, 1 );
		core::Vector p2( 7, 7, 10);
		
		core::Vector q1( 3, 4, 5 );
		core::Vector q2( 8, 3, 1 );
		
		TS_ASSERT_EQUALS( clockwise( p1, p2, q1 ), false );
		TS_ASSERT_EQUALS( clockwise( p1, p2, q2 ), true );
		TS_ASSERT_EQUALS( clockwise( p2, p1, q1 ), true );
		TS_ASSERT_EQUALS( clockwise( p2, p1, q2 ), false );
		
		
	}
	
	// test distance calculation of q from line p1p2 in 2D
	void test_distance() {
		
		TR << "test distance" << std::endl;
		
		core::Vector p1( 1, 6, 5 );
		core::Vector p2( 7, 6, 8 );
		core::Vector p3( 7, 10, 2 );

		core::Vector q1( 5, 9, 2 );
		core::Vector q2( 2, 7, 10 );
		
		TS_ASSERT_EQUALS( distance_from_line2D( p1, p2, q1 ), 3 );
		TS_ASSERT_EQUALS( distance_from_line2D( p2, p3, q2 ), 5 );
		
	}

	// test output vector of points that are in triangle
	void test_points_in_triangle() {
		
		TR << "test points in triangle" << std::endl;
		
		// triangle
		core::Vector p1( 1, 1, 1 );
		core::Vector p2( 7, 0, 10);
		core::Vector p3( 7, 7, 7 );

		// points inside
		core::Vector q1( 3, 2, 8 );
		core::Vector q2( 6, 3, 3 );

		// points outside
		core::Vector q3( 3, 5, 12 );
		core::Vector q4( 9, 4, 1 );

		// triangle points
		std::map< core::Size, core::Vector > coords;
		coords[ 1 ] = p1;
		coords[ 2 ] = p2;
		coords[ 3 ] = p3;

		// inside and outside
		coords[ 4 ] = q1;
		coords[ 5 ] = q2;
		coords[ 6 ] = q3;
		coords[ 7 ] = q4;
		
		// pointlist to test
		utility::vector1< core::Size > pointlist;
		pointlist.push_back( 4 );
		pointlist.push_back( 5 );
		pointlist.push_back( 6 );
		pointlist.push_back( 7 );
		
		// output
		utility::vector1< core::Size > output;
		output.push_back( 4 );
		output.push_back( 5 );
		
		// 1, 2, 3 are point numbers for the triangle
		TS_ASSERT_EQUALS( points_in_triangle( coords, 1, 2, 3, pointlist ), output );
		
	}

	// test enclosing angles
	void test_enclosing_angles() {
		
		TR << "test enclosing angles" << std::endl;
		
		core::Vector p1( 1, 1, 1 );
		core::Vector p2( 5, 5, 5 );
		core::Vector q( 5, 1, 2 );
		
		TS_ASSERT_DELTA( enclosing_angles( p1, p2, q ), 90.0, 0.0001 );
		
	}
	
	// get distances
	void test_get_distances() {
		
		TR << "test get distances" << std::endl;
		
		// points
		core::Vector q1( 1, 1, 1 );
		core::Vector q2( 7, 0, 10);
		core::Vector q3( 7, 7, 7 );
		core::Vector q4( -3, -2, 8 );
		
		core::Vector p1( 2, -2, 3 );
		core::Vector p2( -2, 2, 9 );
		
		// list
		utility::vector1< core::Size > outside;
		outside.push_back( 1 );
		outside.push_back( 2 );
		outside.push_back( 3 );
		outside.push_back( 4 );
		
		// coordinates
		std::map< core::Size, core::Vector > coords;
		coords[ 1 ] = q1;
		coords[ 2 ] = q2;
		coords[ 3 ] = q3;
		coords[ 4 ] = q4;
		
		coords[ 5 ] = p1;
		coords[ 6 ] = p2;
		
		utility::vector1< core::Real > distances = get_distances( coords, 5, 6, outside, false );

		TS_ASSERT_DELTA( distances[1], 1.41421, 0.0001 );
		TS_ASSERT_DELTA( distances[2], 4.94975, 0.0001 );
		TS_ASSERT_DELTA( distances[3], 9.89949, 0.0001 );
		TS_ASSERT_DELTA( distances[4], -1000.0, 0.0001 );
		
	}

	// test find farthest
	void test_find_farthest() {
		
		TR << "test find farthest" << std::endl;
		
		// points
		core::Vector q1( 1, 1, 1 );
		core::Vector q2( 7, 0, 10);
		core::Vector q3( 7, 7, 7 );
		core::Vector q4( 3, 2, 8 );
		core::Vector q5( 6, 3, 3 );
		core::Vector q6( 3, 5, 12 );
		core::Vector q7( 9, 4, 1 );

		core::Vector p1( 2, -2, 3 );
		core::Vector p2( -2, 2, 9 );
		
		// list
		utility::vector1< core::Size > outside;
		outside.push_back( 1 );
		outside.push_back( 2 );
		outside.push_back( 3 );
		outside.push_back( 4 );
		outside.push_back( 5 );
		outside.push_back( 6 );
		outside.push_back( 7 );
		
		// coordinates
		std::map< core::Size, core::Vector > coords;
		coords[ 1 ] = q1;
		coords[ 2 ] = q2;
		coords[ 3 ] = q3;
		coords[ 4 ] = q4;
		coords[ 5 ] = q5;
		coords[ 6 ] = q6;
		coords[ 7 ] = q7;
		coords[ 8 ] = p1;
		coords[ 9 ] = p2;
		
		TS_ASSERT_EQUALS( find_farthest( coords, 8, 9, outside, false ), 3 );

	}
	
	// get angles
	void test_get_angles() {
		
		TR << "test get angles" << std::endl;
		
		// points
		core::Vector q1( 1, 1, 1 );
		core::Vector q2( 7, 0, 10);
		core::Vector q3( 7, 7, 7 );
		core::Vector q4( -3, -2, 8 );
		
		core::Vector p1( 2, -2, 3 );
		core::Vector p2( -2, 2, 9 );
		
		// list
		utility::vector1< core::Size > outside;
		outside.push_back( 1 );
		outside.push_back( 2 );
		outside.push_back( 3 );
		outside.push_back( 4 );
		
		// coordinates
		std::map< core::Size, core::Vector > coords;
		coords[ 1 ] = q1;
		coords[ 2 ] = q2;
		coords[ 3 ] = q3;
		coords[ 4 ] = q4;
		
		coords[ 5 ] = p1;
		coords[ 6 ] = p2;
		
		utility::vector1< core::Real > angles = get_angles( coords, 5, 6, outside, false );
		
		TS_ASSERT_DELTA( angles[1], 53.1301, 0.0001 );
		TS_ASSERT_DELTA( angles[2], 145.6697, 0.0001 );
		TS_ASSERT_DELTA( angles[3], 148.1092, 0.0001 );
		TS_ASSERT_DELTA( angles[4], 9999, 0.0001 );
		
	}

	// find closest
	void test_find_closest() {
		
		TR << "test find closest" << std::endl;
		
		// points
		core::Vector q1( 1, 1, 1 );
		core::Vector q2( 7, 0, 10);
		core::Vector q3( 7, 7, 7 );
		core::Vector q4( 3, 2, 8 );
		core::Vector q5( 6, 3, 3 );
		core::Vector q6( 3, 5, 12 );
		core::Vector q7( 9, 4, 1 );
		
		core::Vector p1( 2, -2, 3 );
		core::Vector p2( -2, 2, 9 );
		
		// list
		utility::vector1< core::Size > outside;
		outside.push_back( 1 );
		outside.push_back( 2 );
		outside.push_back( 3 );
		outside.push_back( 4 );
		outside.push_back( 5 );
		outside.push_back( 6 );
		outside.push_back( 7 );
		
		// coordinates
		std::map< core::Size, core::Vector > coords;
		coords[ 1 ] = q1;
		coords[ 2 ] = q2;
		coords[ 3 ] = q3;
		coords[ 4 ] = q4;
		coords[ 5 ] = q5;
		coords[ 6 ] = q6;
		coords[ 7 ] = q7;
		coords[ 8 ] = p1;
		coords[ 9 ] = p2;
		
		TS_ASSERT_EQUALS( core::membrane::find_closest( coords, 8, 9, outside, false ), 1 );
		
	}
	
};
