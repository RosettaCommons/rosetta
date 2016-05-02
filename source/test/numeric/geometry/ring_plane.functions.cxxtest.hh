// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/numeric/geometry/ring_plane.functions.cxxtest.hh
/// @brief   Test suite for ring-plane geometry functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <core/init_util.hh>

// Unit header
#include <numeric/geometry/ring_plane.functions.hh>

// Numeric headers
#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.hh>
#include <numeric/xyz.functions.hh>

// Utility header
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <basic/Tracer.hh>

// C++ header
#include <math.h>

static THREAD_LOCAL basic::Tracer TR("numeric.geometry.ring_plane.functions.cxxtest");

class RingPlaneFunctionsTests : public CxxTest::TestSuite {
public: // Type definitions ///////////////////////////////////////////////////
	typedef numeric::xyzVector< numeric::Real > Coords;


public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init(); // So that tracers get muted
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	// Confirm that R-squared values are properly reported for points relative to a plane.
	void test_residual_squared_of_points_to_plane()
	{
		using namespace numeric;
		using namespace numeric::geometry;
		using namespace utility;

		TR <<  "Testing that residual_squared_of_points_to_plane() properly calculates R-squared values."  << std::endl;

		vector1< Coords > point( 1 ), segment( 2 ), square( 4 );
		xyzVector< Real > const xy_plane( 0.0, 0.0, 1.0 );

		point[ 1 ] = Coords( 0.0 );  // point at origin is on the xy plane
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( point, xy_plane ), 0.0 );

		point[ 1 ] = Coords( 1.0, 0.0, 0.0 );  // point at x=1...
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( point, xy_plane ), 0.0 );

		point[ 1 ] = Coords( 0.0, 0.0, 1.0 );  // point at z=1 is off the xy plane, but a point is always "parallel"
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( point, xy_plane ), 0.0 );


		segment[ 1 ] = Coords( -1.0, 0.0, 0.0 );  // segment along x axis is on the xy plane
		segment[ 2 ] = Coords( 1.0, 0.0, 0.0 );
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( segment, xy_plane ), 0.0 );

		segment[ 1 ] = Coords( 0.0, 0.0, 0.0 );  // segment between x and z axes
		segment[ 2 ] = Coords( 1.0, 0.0, 1.0 );
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( segment, xy_plane ), 0.5 );

		segment[ 1 ] = Coords( 0.0, 0.0, -1.0 );  // segment along z axis is normal to the xy plane
		segment[ 2 ] = Coords( 0.0, 0.0, 1.0 );
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( segment, xy_plane ), 2.0 );


		square[ 1 ] = Coords( -1.0, 0.0, 0.0 );  // square in the xy plane (like a diamond)
		square[ 2 ] = Coords( 0.0, 1.0, 0.0 );
		square[ 3 ] = Coords( 1.0, 0.0, 0.0 );
		square[ 4 ] = Coords( 0.0, -1.0, 0.0 );
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( square, xy_plane ), 0.0 );

		square[ 1 ] = Coords( -1.0, 0.0, 1.0 );  // square parallel to the xy plane, so R^2 is still 0
		square[ 2 ] = Coords( 0.0, 1.0, 1.0 );
		square[ 3 ] = Coords( 1.0, 0.0, 1.0 );
		square[ 4 ] = Coords( 0.0, -1.0, 1.0 );
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( square, xy_plane ), 0.0 );

		square[ 1 ] = Coords( 0.0, -1.0, 0.0 );  // square in the yz plane bisecting the xy plane
		square[ 2 ] = Coords( 0.0, 0.0, 1.0 );
		square[ 3 ] = Coords( 0.0, 1.0, 0.0 );
		square[ 4 ] = Coords( 0.0, 0.0, -1.0 );
		TS_ASSERT_EQUALS( residual_squared_of_points_to_plane( square, xy_plane ), 2.0 );
	}


	// Confirm that the coordinates of the points of the planar ring system lie in the calculated average plane of that
	// ring system.
	void test_planar_rings_without_checking_for_atom_planarity()
	{
		using namespace numeric;
		using namespace numeric::geometry;
		using namespace numeric::constants::r;
		using namespace utility;

		// Constants
		Real const sin60( sin( 60.0 * pi_over_180 ) );  // = square root of 3 over 2
		Real const sin45over2( sin( 45.0 * pi_over_180 ) / 2 );  // = square root of 2 over 4
		Real const sin45times2( sin( 45.0 * pi_over_180 ) * 2 );  // = square root of 2


		TR <<  "Testing how vector_normal_to_ring_plane_of_best_fit() handles non- and planar rings."  << std::endl;

		vector1< Coords > ring1( 1 ), ring2( 2 ), ring3( 3 ), /*ring4( 4 ), ring5( 5 ),*/ ring6( 6 );

		// 1-membered "ring"
		ring1[ 1 ] = Coords( 0.0 );  // just a point

		TS_ASSERT( are_coplanar( ring1 ) );
		TS_ASSERT_EQUALS( vector_normal_to_ring_plane_of_best_fit( ring1, false ), ZERO_VECTOR );

		// 2-membered "ring"
		ring2[ 1 ] = Coords( 0.0 );
		ring2[ 2 ] = Coords( 1.0 );  // just a line segment

		TS_ASSERT( are_coplanar( ring2 ) );
		TS_ASSERT_EQUALS( vector_normal_to_ring_plane_of_best_fit( ring2, false ), ZERO_VECTOR );

		// 3-membered ring (cyclopropane carbons, coordinates generated by Discovery Studio.)
		ring3[ 1 ] = Coords( -1.114, 0.281, -0.404 );
		ring3[ 2 ] = Coords( -0.120, -0.895, -0.404 );
		ring3[ 3 ] = Coords( 0.398, 0.551, -0.296 );

		TS_ASSERT( are_coplanar( ring3 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring3, vector_normal_to_ring_plane_of_best_fit( ring3, false ) ), 0.0, 0.02 );

		// 6-membered ring (ideal hexagon stretched and tilted out of the xy plane by raising or lowering two vertical
		// sides by 1 unit in the z directions)
		ring6[ 1 ] = Coords( 0.0, 1.0, 0.0 );
		ring6[ 2 ] = Coords( sin60, 0.5, 1.0 );
		ring6[ 3 ] = Coords( sin60, -0.5, 1.0 );
		ring6[ 4 ] = Coords( 0.0, -1.0, 0.0 );
		ring6[ 5 ] = Coords( -sin60, -0.5, -1.0 );
		ring6[ 6 ] = Coords( -sin60, 0.5, -1.0 );

		TS_ASSERT( are_coplanar( ring6 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6, false ) ), 0.0, 0.02 );


		// 6-membered ring in a chair conformation
		ring6[ 1 ] = Coords( 0.0, 1.0, sin45over2 );
		ring6[ 2 ] = Coords( sin60, 0.5, 0.0 );
		ring6[ 3 ] = Coords( sin60, -0.5, sin45over2 );
		ring6[ 4 ] = Coords( 0.0, -1.0, 0.0 );
		ring6[ 5 ] = Coords( -sin60, -0.5, sin45over2 );
		ring6[ 6 ] = Coords( -sin60, 0.5, 0.0 );
		TS_ASSERT( ! are_coplanar( ring6 ) );

		// Since points 2, 4, and 6 lie in the xy plane, the average plane parallel to xy and bisecting this ring
		// system must have the origin of its normal vector shifted by a distance of half that of a tetrahedral
		// inradius in the z direction, or 1/root 32 units.  The z component of the displacement vectors of each point
		// in the ring system from the center of mass would be plus or minus 1/root 32.  So R-squared will be 6 times
		// 1/32 = 3/16 = 0.1875.
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6, false ) ), 0.1875, 0.02 );

		// the same 6-membered ring in a chair conformation rotated randomly using DiscoveryStudio
		ring6[ 1 ] = Coords( -0.9812, -0.0692, 0.2524 );
		ring6[ 2 ] = Coords( -0.2712, -0.7819, 0.5885 );
		ring6[ 3 ] = Coords( 0.5064, -0.8709, -0.1274 );
		ring6[ 4 ] = Coords( 0.9812, 0.0692, -0.2524 );
		ring6[ 5 ] = Coords( 0.2712, 0.7819, -0.5885 );
		ring6[ 6 ] = Coords( -0.5064, 0.8709, 0.1274 );

		TS_ASSERT( ! are_coplanar( ring6 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6, false ) ), 0.1875, 0.02 );


		// 6-membered ring in a boat conformation
		ring6[ 1 ] = Coords( 0.0, sin45over2, 0.5 );
		ring6[ 2 ] = Coords( sin60, 0.0, 0.0 );
		ring6[ 3 ] = Coords( sin60, -3 * sin45over2, 0 );
		ring6[ 4 ] = Coords( 0.0, -sin45times2, 0.5 );
		ring6[ 5 ] = Coords( -sin60, -3 * sin45over2, 0 );
		ring6[ 6 ] = Coords( -sin60, 0.0, 0.0 );

		// Since points 2, 3, 5, and 6 lie in the xy plane, the average plane parallel to xy and bisecting this ring
		// system must have the origin of its normal vector shifted by a distance of 1/6 units in the z direction.
		// The z component of the displacement vectors of points 1 and 4 from the center of mass would be 1/2 − 1/6 =
		// 1/3, and the z component of the displacement vectors of the other points would be 0 − 1/6 = 1/6.  So R2 will
		// be 2 times 1/9 + 4 times 1/36 = 1/3 = 0.3333.

		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6, false ) ), 0.3333, 0.002 );

		TS_ASSERT( ! are_coplanar( ring6 ) );


		// the same 6-membered ring in a boat conformation rotated and translated randomly using DiscoveryStudio
		ring6[ 1 ] = Coords( 1.0313, 2.0483, 0.2160 );
		ring6[ 2 ] = Coords( 0.9235, 1.0355, 0.5120 );
		ring6[ 3 ] = Coords( 1.3354, 0.4033, -0.2335 );
		ring6[ 4 ] = Coords( 1.7179, 0.9947, -1.0265 );
		ring6[ 5 ] = Coords( 2.4725, 1.6403, -0.6540 );
		ring6[ 6 ] = Coords( 2.0605, 2.2724, 0.0915 );

		TS_ASSERT ( ! are_coplanar( ring6 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6, false ) ), 0.3333, 0.05 );

	}

	// Confirm that the coordinates of the points of the planar ring system lie in the calculated average plane of that
	// ring system.
	void test_planar_rings()
	{
		using namespace numeric;
		using namespace numeric::geometry;
		using namespace utility;
		// Constant
		Real const sin60( sin( 60.0 * constants::r::pi_over_180 ) );  // = square root of 3 over 2


		TR <<  "Testing how vector_normal_to_ring_plane_of_best_fit() handles non- and planar rings."  << std::endl;

		vector1< Coords > ring1( 1 ), ring2( 2 ), ring3( 3 ), /*ring4( 4 ), ring5( 5 ),*/ ring6( 6 );

		// 1-membered "ring"
		ring1[ 1 ] = Coords( 0.0 );  // just a point

		TS_ASSERT( are_coplanar( ring1 ) );
		TS_ASSERT_EQUALS( vector_normal_to_ring_plane_of_best_fit( ring1 ), ZERO_VECTOR );

		// 2-membered "ring"
		ring2[ 1 ] = Coords( 0.0 );
		ring2[ 2 ] = Coords( 1.0 );  // just a line segment

		TS_ASSERT( are_coplanar( ring2 ) );
		TS_ASSERT_EQUALS( vector_normal_to_ring_plane_of_best_fit( ring2 ), ZERO_VECTOR );

		// 3-membered ring (cyclopropane carbons, coordinates generated by Discovery Studio.)
		ring3[ 1 ] = Coords( -1.114, 0.281, -0.404 );
		ring3[ 2 ] = Coords( -0.120, -0.895, -0.404 );
		ring3[ 3 ] = Coords( 0.398, 0.551, -0.296 );

		TS_ASSERT( are_coplanar( ring3 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring3, vector_normal_to_ring_plane_of_best_fit( ring3 ) ), 0.0, 0.02 );

		// 6-membered ring (ideal hexagon stretched and tilted out of the xy plane by raising or lowering two vertical
		// sides by 1 unit in the z directions)
		ring6[ 1 ] = Coords( 0.0, 1.0, 0.0 );
		ring6[ 2 ] = Coords( sin60, 0.5, 1.0 );
		ring6[ 3 ] = Coords( sin60, -0.5, 1.0 );
		ring6[ 4 ] = Coords( 0.0, -1.0, 0.0 );
		ring6[ 5 ] = Coords( -sin60, -0.5, -1.0 );
		ring6[ 6 ] = Coords( -sin60, 0.5, -1.0 );

		TS_ASSERT( are_coplanar( ring6 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6 ) ), 0.0, 0.02 );
	}

	// Confirm that non-planar rings have the expected residuals from their calculated average planes.
	void test_non_planar_rings()
	{
		using namespace numeric;
		using namespace numeric::geometry;
		using namespace numeric::constants::r;
		using namespace utility;

		// Constants
		Real const sin60( sin( 60.0 * pi_over_180 ) );  // = square root of 3 over 2
		Real const sin45over2( sin( 45.0 * pi_over_180 ) / 2 );  // = square root of 2 over 4
		Real const sin45times2( sin( 45.0 * pi_over_180 ) * 2 );  // = square root of 2


		TR <<  "Testing how vector_normal_to_ring_plane_of_best_fit() handles non-planar rings."  << std::endl;

		vector1< Coords > ring6( 6 );

		// 6-membered ring in a chair conformation
		ring6[ 1 ] = Coords( 0.0, 1.0, sin45over2 );
		ring6[ 2 ] = Coords( sin60, 0.5, 0.0 );
		ring6[ 3 ] = Coords( sin60, -0.5, sin45over2 );
		ring6[ 4 ] = Coords( 0.0, -1.0, 0.0 );
		ring6[ 5 ] = Coords( -sin60, -0.5, sin45over2 );
		ring6[ 6 ] = Coords( -sin60, 0.5, 0.0 );
		TS_ASSERT( ! are_coplanar( ring6 ) );

		// Since points 2, 4, and 6 lie in the xy plane, the average plane parallel to xy and bisecting this ring
		// system must have the origin of its normal vector shifted by a distance of half that of a tetrahedral
		// inradius in the z direction, or 1/root 32 units.  The z component of the displacement vectors of each point
		// in the ring system from the center of mass would be plus or minus 1/root 32.  So R-squared will be 6 times
		// 1/32 = 3/16 = 0.1875.
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6 ) ), 0.1875, 0.02 );

		// the same 6-membered ring in a chair conformation rotated randomly using DiscoveryStudio
		ring6[ 1 ] = Coords( -0.9812, -0.0692, 0.2524 );
		ring6[ 2 ] = Coords( -0.2712, -0.7819, 0.5885 );
		ring6[ 3 ] = Coords( 0.5064, -0.8709, -0.1274 );
		ring6[ 4 ] = Coords( 0.9812, 0.0692, -0.2524 );
		ring6[ 5 ] = Coords( 0.2712, 0.7819, -0.5885 );
		ring6[ 6 ] = Coords( -0.5064, 0.8709, 0.1274 );

		TS_ASSERT( ! are_coplanar( ring6 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6 ) ), 0.1875, 0.02 );


		// 6-membered ring in a boat conformation
		ring6[ 1 ] = Coords( 0.0, sin45over2, 0.5 );
		ring6[ 2 ] = Coords( sin60, 0.0, 0.0 );
		ring6[ 3 ] = Coords( sin60, -3 * sin45over2, 0 );
		ring6[ 4 ] = Coords( 0.0, -sin45times2, 0.5 );
		ring6[ 5 ] = Coords( -sin60, -3 * sin45over2, 0 );
		ring6[ 6 ] = Coords( -sin60, 0.0, 0.0 );

		// Since points 2, 3, 5, and 6 lie in the xy plane, the average plane parallel to xy and bisecting this ring
		// system must have the origin of its normal vector shifted by a distance of 1/6 units in the z direction.
		// The z component of the displacement vectors of points 1 and 4 from the center of mass would be 1/2 − 1/6 =
		// 1/3, and the z component of the displacement vectors of the other points would be 0 − 1/6 = 1/6.  So R2 will
		// be 2 times 1/9 + 4 times 1/36 = 1/3 = 0.3333.

		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6 ) ), 0.3333, 0.002 );

		TS_ASSERT( ! are_coplanar( ring6 ) );


		// the same 6-membered ring in a boat conformation rotated and translated randomly using DiscoveryStudio
		ring6[ 1 ] = Coords( 1.0313, 2.0483, 0.2160 );
		ring6[ 2 ] = Coords( 0.9235, 1.0355, 0.5120 );
		ring6[ 3 ] = Coords( 1.3354, 0.4033, -0.2335 );
		ring6[ 4 ] = Coords( 1.7179, 0.9947, -1.0265 );
		ring6[ 5 ] = Coords( 2.4725, 1.6403, -0.6540 );
		ring6[ 6 ] = Coords( 2.0605, 2.2724, 0.0915 );

		TS_ASSERT ( ! are_coplanar( ring6 ) );
		TS_ASSERT_DELTA( residual_squared_of_points_to_plane(
				ring6, vector_normal_to_ring_plane_of_best_fit( ring6 ) ), 0.3333, 0.05 );
	}
};  // class RingPlaneFunctionsTests
