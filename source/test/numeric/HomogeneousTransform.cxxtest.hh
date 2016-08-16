// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/HomogeneousTransform.cxxtest.hh
/// @brief  test suite for numeric::HomogeneousTransform.cxxtest.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyzVector.io.hh>

typedef numeric::HomogeneousTransform< double > HTD;

using namespace numeric;

// --------------- Test Class --------------- //

class HomogenousTransformTests : public CxxTest::TestSuite {

	public:
		typedef xyzVector< double > Vector;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	/// @brief Test identity transform constructor
	void test_HomogeneousTransform_default_constructor() {
		HTD ht;
		TS_ASSERT_DELTA( ht.xx(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( ht.xy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.yx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.yy(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( ht.yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.zz(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( ht.px(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.py(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( ht.pz(), 0.0, 1e-14 );
	}

	void test_HomogeneousTransform_matrix_constructor() {
		xyzMatrix< double > zrotation = z_rotation_matrix_degrees( 30 );

		xyzVector< double > point( 1.5, 2.25, -3.125 );

		double const sin30 = sin( conversions::radians( 30 ) );
		double const cos30 = cos( conversions::radians( 30 ) );

		/// Double check that z_rotation_matrix is working properly;
		TS_ASSERT_DELTA( zrotation.xx(), cos30, 1e-14 );
		TS_ASSERT_DELTA( zrotation.yx(), -sin30, 1e-14 );
		TS_ASSERT_DELTA( zrotation.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( zrotation.xy(), sin30, 1e-14 );
		TS_ASSERT_DELTA( zrotation.yy(), cos30, 1e-14 );
		TS_ASSERT_DELTA( zrotation.zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( zrotation.xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( zrotation.yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( zrotation.zz(), 1.0, 1e-14 );


		HTD z_ht( zrotation, point );
		TS_ASSERT_DELTA( z_ht.xx(), cos30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.xy(), -sin30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.yx(), sin30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.yy(), cos30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.zz(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.px(), 1.5, 1e-14 );
		TS_ASSERT_DELTA( z_ht.py(), 2.25, 1e-14 );
		TS_ASSERT_DELTA( z_ht.pz(), -3.125, 1e-14 );

		// Check that rotation_matrix accessor generates correct matrix
		TS_ASSERT_DELTA( z_ht.rotation_matrix().xx(), cos30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().yx(), -sin30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().xy(), sin30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().yy(), cos30, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( z_ht.rotation_matrix().zz(), 1.0, 1e-14 );

		xyzMatrix< double > full_rotation = (
				z_rotation_matrix_degrees( 30 ) *
				y_rotation_matrix_degrees( 20 ) *
				x_rotation_matrix_degrees( 10 ));

		HTD full_ht( zrotation, point );

		// Check that rotation_matrix accessor generates correct matrix
		TS_ASSERT_DELTA( full_ht.rotation_matrix().xx(), full_rotation.xx(), 1e-14 );
		TS_ASSERT_DELTA( full_ht.rotation_matrix().xy(), full_rotation.xy(), 1e-14 );
		TS_ASSERT_DELTA( full_ht.rotation_matrix().xz(), full_rotation.xz(), 1e-14 );

		TS_ASSERT_DELTA( full_ht.rotation_matrix().yx(), full_rotation.yx(), 1e-14 );
		TS_ASSERT_DELTA( full_ht.rotation_matrix().yy(), full_rotation.yy(), 1e-14 );
		TS_ASSERT_DELTA( full_ht.rotation_matrix().yz(), full_rotation.yz(), 1e-14 );
		
		TS_ASSERT_DELTA( full_ht.rotation_matrix().zx(), full_rotation.zx(), 1e-14 );
		TS_ASSERT_DELTA( full_ht.rotation_matrix().zy(), full_rotation.zy(), 1e-14 );
		TS_ASSERT_DELTA( full_ht.rotation_matrix().zz(), full_rotation.zz(), 1e-14 );
	}

	void test_HomogenousTransform_multiplication() {

		/// 1. Walk along the zaxis a distance of 1.5.
		/// The new coordinate frame will still have the identity transform,
		/// but it's coordinate will be at (0,0,1.5).
		/// Right multiply the series of motions to create the forward kinematics
		/// calculation of where we end up.
		HTD ht;
		HTD walk_1p5z;
		walk_1p5z.set_transform( xyzVector< double >( 0.0, 0.0, 1.5 ) );
		HTD newloc = ht * walk_1p5z;
		TS_ASSERT_DELTA( newloc.xx(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.xy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.yx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.yy(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.zz(), 1.0, 1e-14 );

		TS_ASSERT_DELTA( newloc.px(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.py(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc.pz(), 1.5, 1e-14 );


		//// 2. Perform a negative rotation about the xaxis to mimic a bond angle
		//// and then walk from there another 1.5 A.
		HTD twist_x70p5;
		/// tetrahedral geometry (180 - 109.5 = 70.5; negative rotation about x axis to walk into the first quadrant of the yz plane )
		twist_x70p5.set_xaxis_rotation_deg( -70.5 );
		HTD newloc2 =  newloc * twist_x70p5 * walk_1p5z;

		double const sin19p5 = sin( conversions::radians( 19.5 ) );
		double const cos19p5 = cos( conversions::radians( 19.5 ) );

		TS_ASSERT_DELTA( newloc2.px(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc2.py(), 1.5 * cos19p5, 1e-14 );
		TS_ASSERT_DELTA( newloc2.pz(), 1.5 + 1.5 * sin19p5, 1e-14 );

		/// test x axis still entirely in x plane
		TS_ASSERT_DELTA( newloc2.xx(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( newloc2.xy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc2.xz(), 0.0, 1e-14 );

		/// test y axis now in yz plane:
		TS_ASSERT_DELTA( newloc2.yx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc2.yy(), sin19p5, 1e-14 );
		TS_ASSERT_DELTA( newloc2.yz(), -cos19p5, 1e-14 );

		/// test zaxis is also in yz plane
		TS_ASSERT_DELTA( newloc2.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc2.zy(), cos19p5, 1e-14 );
		TS_ASSERT_DELTA( newloc2.zz(), sin19p5, 1e-14 );

		HTD twist_walk_premultiply = twist_x70p5 * walk_1p5z;
		HTD newloc2prime =  newloc * twist_walk_premultiply;

		TS_ASSERT_DELTA( newloc2prime.px(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( newloc2prime.py(), 1.5 * cos19p5, 1e-14 );
		TS_ASSERT_DELTA( newloc2prime.pz(), 1.5 + 1.5 * sin19p5, 1e-14 );

		/// 3. Perform a rotation about the Z axis, mimicing a chi dihedral rotation
		/// before walking to the third point.
		HTD twist_z60;
		twist_z60.set_zaxis_rotation_deg( 60.0 );
		HTD twist_twist_walk_premultiply =  twist_z60 * twist_walk_premultiply;
		HTD newloc3 = newloc2 * twist_twist_walk_premultiply;

		double const sin60 = sin( conversions::radians( 60.0 ) );
		double const cos60 = cos( conversions::radians( 60.0 ) );

		TS_ASSERT_DELTA( twist_z60.xx(), cos60, 1e-14 );
		TS_ASSERT_DELTA( twist_z60.xy(), sin60, 1e-14 );
		TS_ASSERT_DELTA( twist_z60.xz(),     0, 1e-14 );

		TS_ASSERT_DELTA( twist_z60.yx(), -sin60, 1e-14 );
		TS_ASSERT_DELTA( twist_z60.yy(),  cos60, 1e-14 );
		TS_ASSERT_DELTA( twist_z60.yz(),      0, 1e-14 );

		TS_ASSERT_DELTA( twist_z60.zx(), 0, 1e-14 );
		TS_ASSERT_DELTA( twist_z60.zy(), 0, 1e-14 );
		TS_ASSERT_DELTA( twist_z60.zz(), 1, 1e-14 );


		/// test x axis
		TS_ASSERT_DELTA( newloc3.xx(),          cos60, 1e-14 );
		TS_ASSERT_DELTA( newloc3.xy(),  sin60*sin19p5, 1e-14 );
		TS_ASSERT_DELTA( newloc3.xz(), -sin60*cos19p5, 1e-14 );

		/// test y axis
		TS_ASSERT_DELTA( newloc3.yx(), -sin19p5*sin60,                           1e-14 );
		TS_ASSERT_DELTA( newloc3.yy(), -cos19p5*cos19p5 + sin19p5*cos60*sin19p5, 1e-14 );
		TS_ASSERT_DELTA( newloc3.yz(), -sin19p5*cos19p5 - sin19p5*cos60*cos19p5, 1e-14 );

		/// test z axis
		TS_ASSERT_DELTA( newloc3.zx(), -cos19p5 * sin60,                               1e-14 );
		TS_ASSERT_DELTA( newloc3.zy(),  cos19p5 * cos60 * sin19p5 + sin19p5 * cos19p5, 1e-14 );
		TS_ASSERT_DELTA( newloc3.zz(), -cos19p5 * cos60 * cos19p5 + sin19p5 * sin19p5, 1e-14 );

		/// test coordinate
		TS_ASSERT_DELTA( newloc3.px(), newloc2.px() + 1.5 * (-cos19p5 * sin60),                               1e-14 );
		TS_ASSERT_DELTA( newloc3.py(), newloc2.py() + 1.5 * ( cos19p5 * cos60 * sin19p5 + sin19p5 * cos19p5), 1e-14 );
		TS_ASSERT_DELTA( newloc3.pz(), newloc2.pz() + 1.5 * (-cos19p5 * cos60 * cos19p5 + sin19p5 * sin19p5), 1e-14 );

	}

	void test_HomogenousTransform_vector_multiplication() {

		/// calculate the location of a point in the global coordinate frame
		/// given its location in the coordinate frame constructed from newloc2

		/// 1. Walk along the zaxis a distance of 1.5.
		/// The new coordinate frame will still have the identity transform,
		/// but it's coordinate will be at (0,0,1.5).
		HTD walk_1p5z;
		walk_1p5z.set_transform( xyzVector< double >( 0.0, 0.0, 1.5 ) );

		//// 2. Perform a negative rotation about the xaxis to mimic a bond angle
		//// and then walk from there another 1.5 A.
		HTD twist_x70p5;
		/// tetrahedral geometry (180 - 109.5 = 70.5; negative rotation about x axis to walk into the first quadrant of the yz plane )
		twist_x70p5.set_xaxis_rotation_deg( -70.5 );
		HTD newloc2 = /*start at the origin, and */ walk_1p5z * twist_x70p5 * walk_1p5z;

		HTD twist_walk_premultiply = twist_x70p5 * walk_1p5z;

		/// 3. Perform a rotation about the Z axis, mimicing a chi dihedral rotation
		/// before walking to the third point.
		HTD twist_z60;
		twist_z60.set_zaxis_rotation_deg( 60.0 );
		HTD twist_twist_walk_premultiply =  twist_z60 * twist_walk_premultiply;

		double const sin19p5 = sin( conversions::radians( 19.5 ) );
		double const cos19p5 = cos( conversions::radians( 19.5 ) );

		double const sin60 = sin( conversions::radians( 60.0 ) );
		double const cos60 = cos( conversions::radians( 60.0 ) );

		/// If we have computed the location of a point in a particular coordinate frame,
		/// then we may construct that point's location given the coordinate frame by doing
		/// a matrix-vector multiply.
		xyzVector< double > point_in_p2 = twist_twist_walk_premultiply.point();

		//std::cout << "point in p2: " << point_in_p2.x() << " " << point_in_p2.y() << " " << point_in_p2.z() << std::endl;

		xyzVector< double > point_in_global_coords = newloc2 *  point_in_p2;

		/// test coordinate
		TS_ASSERT_DELTA( point_in_global_coords.x(), newloc2.px() + 1.5 * (-cos19p5 * sin60),                               1e-14 );
		TS_ASSERT_DELTA( point_in_global_coords.y(), newloc2.py() + 1.5 * ( cos19p5 * cos60 * sin19p5 + sin19p5 * cos19p5), 1e-14 );
		TS_ASSERT_DELTA( point_in_global_coords.z(), newloc2.pz() + 1.5 * (-cos19p5 * cos60 * cos19p5 + sin19p5 * sin19p5), 1e-14 );

	}

	void test_euler_angle_measurement() {

		HTD rot1, rot2, rot3, rot4, rot5;

		rot1.set_xaxis_rotation_deg( 30 );
		rot2.set_yaxis_rotation_deg( 26 );
		rot3.set_zaxis_rotation_deg( 40 );
		rot4.set_yaxis_rotation_deg( 33 );
		rot5.set_zaxis_rotation_deg( 61 );

		HTD frame = rot1 * rot2 * rot3 * rot4 * rot5;  /// define an arbitrary frame
		//HTD frame = rot5;
		xyzVector< double > euler = frame.euler_angles_rad();

		//euler( 1 ) = 0;

		xyzVector< double > euler_deg = euler;
		euler_deg *= numeric::constants::d::radians_to_degrees;

		//std::cout << "Euler measured as: " << euler_deg( 1 ) << " " << euler_deg( 2 ) << " " << euler_deg( 3 ) << std::endl;


		HTD frame2;
		frame2.from_euler_angles_rad( euler );


		TS_ASSERT_DELTA( frame.xx(), frame2.xx(), 1e-14 );
		TS_ASSERT_DELTA( frame.xy(), frame2.xy(), 1e-14 );
		TS_ASSERT_DELTA( frame.xz(), frame2.xz(), 1e-14 );

		TS_ASSERT_DELTA( frame.yx(), frame2.yx(), 1e-14 );
		TS_ASSERT_DELTA( frame.yy(), frame2.yy(), 1e-14 );
		TS_ASSERT_DELTA( frame.yz(), frame2.yz(), 1e-14 );

		TS_ASSERT_DELTA( frame.zx(), frame2.zx(), 1e-14 );
		TS_ASSERT_DELTA( frame.zy(), frame2.zy(), 1e-14 );
		TS_ASSERT_DELTA( frame.zz(), frame2.zz(), 1e-14 );


		xyzVector< double > euler2( 40, 61, 30 );
		xyzVector< double > euler2_rad( euler2 );
		euler2_rad *= numeric::constants::d::degrees_to_radians;

		HTD frame3;
		frame3.from_euler_angles_rad( euler2_rad );

		frame = rot3 * rot1 * rot5;

		TS_ASSERT_DELTA( frame.xx(), frame3.xx(), 1e-14 );
		TS_ASSERT_DELTA( frame.xy(), frame3.xy(), 1e-14 );
		TS_ASSERT_DELTA( frame.xz(), frame3.xz(), 1e-14 );

		TS_ASSERT_DELTA( frame.yx(), frame3.yx(), 1e-14 );
		TS_ASSERT_DELTA( frame.yy(), frame3.yy(), 1e-14 );
		TS_ASSERT_DELTA( frame.yz(), frame3.yz(), 1e-14 );

		TS_ASSERT_DELTA( frame.zx(), frame3.zx(), 1e-14 );
		TS_ASSERT_DELTA( frame.zy(), frame3.zy(), 1e-14 );
		TS_ASSERT_DELTA( frame.zz(), frame3.zz(), 1e-14 );

		frame = rot3;
		xyzVector< double > euler3 = frame.euler_angles_deg();
		//std::cout << "frame: xx_ : " << frame.xx() << " " << numeric::constants::d::radians_to_degrees * std::acos( frame.xx() ) << std::endl;

		//std::cout << "euler3: " << euler3( 1 ) << " " << euler3( 2 ) << " " << euler3( 3 ) << std::endl;
		HTD frame4;
		frame4.from_euler_angles_deg( euler3 );
		TS_ASSERT_DELTA( frame.xx(), frame4.xx(), 1e-14 );
		TS_ASSERT_DELTA( frame.xy(), frame4.xy(), 1e-14 );
		TS_ASSERT_DELTA( frame.xz(), frame4.xz(), 1e-14 );

		TS_ASSERT_DELTA( frame.yx(), frame4.yx(), 1e-14 );
		TS_ASSERT_DELTA( frame.yy(), frame4.yy(), 1e-14 );
		TS_ASSERT_DELTA( frame.yz(), frame4.yz(), 1e-14 );

		TS_ASSERT_DELTA( frame.zx(), frame4.zx(), 1e-14 );
		TS_ASSERT_DELTA( frame.zy(), frame4.zy(), 1e-14 );
		TS_ASSERT_DELTA( frame.zz(), frame4.zz(), 1e-14 );

//-0.240611 0.94392 -0.226098 -4.55621
//0.968712 0.218927 -0.116912 3.92414
//-0.0608565 -0.247154 -0.967063 -8.92507

	//  -0.240611328009 0.943920428956 -0.226098236687 -4.55621442642
	//  0.968711862901 0.218926806763 -0.11691184694 3.92414326226
	//  -0.0608565157406 -0.247154358811 -0.967063186877 -8.92507068448

		xyzVector< double > xaxis( -0.240611328009, 0.968711862901, -0.0608565157406 );
		xyzVector< double > yaxis( 0.943920428956, 0.218926806763, -0.247154358811 );
		xyzVector< double > zaxis( -0.226098236687, -0.11691184694, -0.967063186877 );
		xyzVector< double > point( -4.55621, 3.92414, -8.92507 );

		HTD frame5( xaxis, yaxis, zaxis, point );
		xyzVector< double > euler5 = frame5.euler_angles_deg();

		HTD frame5prime;
		frame5prime.from_euler_angles_deg( euler5 );
		frame5prime.set_point( point );

		TS_ASSERT_DELTA( frame5.xx(), frame5prime.xx(), 1e-11 );
		TS_ASSERT_DELTA( frame5.xy(), frame5prime.xy(), 1e-11 );
		TS_ASSERT_DELTA( frame5.xz(), frame5prime.xz(), 1e-11 );

		TS_ASSERT_DELTA( frame5.yx(), frame5prime.yx(), 1e-11 );
		TS_ASSERT_DELTA( frame5.yy(), frame5prime.yy(), 1e-11 );
		TS_ASSERT_DELTA( frame5.yz(), frame5prime.yz(), 1e-11 );

		TS_ASSERT_DELTA( frame5.zx(), frame5prime.zx(), 1e-11 );
		TS_ASSERT_DELTA( frame5.zy(), frame5prime.zy(), 1e-11 );
		TS_ASSERT_DELTA( frame5.zz(), frame5prime.zz(), 1e-11 );

	}

	void test_euler_angle_measurement2() {

		HTD rot1, rot2, rot3, rot4;//, rot5;

		rot1.set_zaxis_rotation_deg( 30 );
		rot2.set_zaxis_rotation_deg( 250 );
		rot3.set_xaxis_rotation_deg( 178 );
		rot4.set_xaxis_rotation_deg( 182 );
		//rot5.set_xaxis_rotation_deg( 50 );

		HTD frame1 = rot1 * rot3 * rot2;
		HTD frame2 = rot1 * rot4 * rot2;
		//HTD frame3 = rot1 * rot5 * rot2;
		//HTD frame = rot5;

		xyzVector< double > euler1 = frame1.euler_angles_rad();
		xyzVector< double > euler2 = frame2.euler_angles_rad();
		//xyzVector< double > euler3 = frame3.euler_angles_rad();


		xyzVector< double > euler_deg1( euler1 ), euler_deg2( euler2 );//, euler_deg3( euler3 );
		euler_deg1 *= numeric::constants::d::radians_to_degrees;
		euler_deg2 *= numeric::constants::d::radians_to_degrees;
		//euler_deg3 *= numeric::constants::d::radians_to_degrees;

		//std::cout << "Euler1 measured as: " << euler_deg1( 1 ) << " " << euler_deg1( 2 ) << " " << euler_deg1( 3 ) << std::endl;
		//std::cout << "Euler2 measured as: " << euler_deg2( 1 ) << " " << euler_deg2( 2 ) << " " << euler_deg2( 3 ) << std::endl;
		//std::cout << "Euler measured as: " << euler_deg3( 1 ) << " " << euler_deg3( 2 ) << " " << euler_deg3( 3 ) << std::endl;
		//std::cout << "Flipped: " << 180 + euler_deg2( 1 ) << " " << 180 + euler_deg2( 2 ) << std::endl;

		/// When phi and psi are "flipped" as z moves from < 180 to > 180,
		/// they are move by 180 degrees
		TS_ASSERT_DELTA( 180 + euler_deg2( 1 ), euler_deg1( 1 ), 1e-6 );
		TS_ASSERT_DELTA( 180 + euler_deg2( 2 ), 360 + euler_deg1( 2 ), 1e-6 );

	}


	void test_ht_inverse() {
		HTD ht1, ht2, ht3;
		ht1.set_zaxis_rotation_deg( 33 );
		ht2.set_yaxis_rotation_deg( 152 );
		ht3.set_xaxis_rotation_deg( 10 );
		HTD frame1 = ht1;
		frame1.walk_along_z( 1.5 );
		HTD frame2 = frame1 * ht2;
		frame2.walk_along_z( 1.3 );
		HTD frame3 = frame2 * ht3;
		frame3.walk_along_z( 1.7 );

		HTD invframe1 = frame1.inverse();
		HTD invframe2 = frame2.inverse();
		HTD invframe3 = frame3.inverse();

		HTD identity1 = frame1 * invframe1;
		HTD identity2 = frame2 * invframe2;
		HTD identity3 = frame3 * invframe3;

		TS_ASSERT_DELTA( identity1.xx(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.xy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.yx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.yy(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.zz(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.px(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.py(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity1.pz(), 0.0, 1e-14 );

		TS_ASSERT_DELTA( identity2.xx(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.xy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.yx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.yy(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.zz(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.px(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.py(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity2.pz(), 0.0, 1e-14 );

		TS_ASSERT_DELTA( identity3.xx(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.xy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.xz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.yx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.yy(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.yz(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.zx(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.zy(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.zz(), 1.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.px(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.py(), 0.0, 1e-14 );
		TS_ASSERT_DELTA( identity3.pz(), 0.0, 1e-14 );

	}

	void test_z_rotation_bug(){
		HTD ht1,ht2;
		ht1.set_zaxis_rotation_deg( 89.0);
		ht2.set_zaxis_rotation_deg(271.0);
		xyzVector< double > euler1 = ht1.euler_angles_deg();
		xyzVector< double > euler2 = ht2.euler_angles_deg();
		TS_ASSERT( euler1.distance(euler2) > 1e-6 );
	}

	/*void sini_match_frame() {
		Vector p1( 5.912,   4.303,   6.525 );
		Vector p2( 7.049,   4.868,   7.060 );
		Vector p3( 7.513,   4.486,   8.322 );

		HTD frame( p1, p2, p3 );
		std::cout << "Global frame" << std::endl;
		std::cout << "  " << frame.xx() << " " << frame.yx() << " " << frame.zx() << " " << frame.px() << std::endl;
		std::cout << "  " << frame.xy() << " " << frame.yy() << " " << frame.zy() << " " << frame.py() << std::endl;
		std::cout << "  " << frame.xz() << " " << frame.yz() << " " << frame.zz() << " " << frame.pz() << std::endl;

		Vector euler = frame.euler_angles_deg();
		std::cout << "Euler angles:";
		for ( int ii = 1; ii <= 3; ++ii ) std::cout << euler( ii ) << " "; std::cout << std::endl;
	}*/

};


