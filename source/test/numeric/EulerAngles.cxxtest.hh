// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/EulerAngles.cxxtest.hh
/// @brief  test suite for numeric::EulerAngles
/// @author Alex Ford (fordas@uw.edu)
// Test headers


#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/EulerAngles.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

using namespace numeric;

class EulerAnglesTests : public CxxTest::TestSuite {

#define TS_ASSERT_MATRIX_EQUALS(m1, m2) \
{\
			TS_ASSERT_DELTA( m1.xx(), m2.xx(), 1e-14 );\
			TS_ASSERT_DELTA( m1.xy(), m2.xy(), 1e-14 );\
			TS_ASSERT_DELTA( m1.xz(), m2.xz(), 1e-14 );\
			TS_ASSERT_DELTA( m1.yx(), m2.yx(), 1e-14 );\
			TS_ASSERT_DELTA( m1.yy(), m2.yy(), 1e-14 );\
			TS_ASSERT_DELTA( m1.yz(), m2.yz(), 1e-14 );\
			TS_ASSERT_DELTA( m1.zx(), m2.zx(), 1e-14 );\
			TS_ASSERT_DELTA( m1.zy(), m2.zy(), 1e-14 );\
			TS_ASSERT_DELTA( m1.zz(), m2.zz(), 1e-14 );\
}

#define TS_ASSERT_EULER_EQUALS(e1, e2) \
{\
	TS_ASSERT_DELTA( e1.phi(),   e2.phi(), 1e-14 );\
	TS_ASSERT_DELTA( e1.psi(),   e2.psi(), 1e-14 );\
	TS_ASSERT_DELTA( e1.theta(), e2.theta(), 1e-14 );\
}

	public:
		typedef xyzVector< double > Vector;
		typedef xyzMatrix< double > Matrix;
		typedef numeric::EulerAngles< double > Angles;

		void test_rotation_matrix_conversion()
		{
			Matrix rot1, rot2, rot3, rot4, rot5;

			rot1 = x_rotation_matrix_degrees( 30. );
			rot2 = y_rotation_matrix_degrees( 26. );
			rot3 = z_rotation_matrix_degrees( 40. );
			rot4 = y_rotation_matrix_degrees( 33. );
			rot5 = z_rotation_matrix_degrees( 61. );

			// Test round-trip conversion from matrix
			Matrix frame =  rot1 * rot2 * rot3 * rot4 * rot5;

			Angles frame_angles(frame);

			Matrix frame_from_angles = frame_angles.to_rotation_matrix();
			TS_ASSERT_MATRIX_EQUALS(frame, frame_from_angles);

			// Test Z-axis rotation
			Matrix z_rotation_frame = rot3;
			Angles z_rotation_euler(0, 0, 0);
			z_rotation_euler.phi_degrees( 40. );

			Matrix z_rotation_frame_from_angles = z_rotation_euler.to_rotation_matrix();
			Angles z_rotation_euler_from_frame(z_rotation_frame);

			TS_ASSERT_EULER_EQUALS(z_rotation_euler, z_rotation_euler_from_frame);
			TS_ASSERT_MATRIX_EQUALS(z_rotation_frame, z_rotation_frame_from_angles);


			// Test init from specific angles
			Angles euler2;
			euler2.phi_degrees( 40. );
			euler2.psi_degrees( 61. );
			euler2.theta_degrees( 30. );

			Matrix frame2 = rot3 * rot1 * rot5;
			Matrix frame2_from_angles = euler2.to_rotation_matrix();

			Angles euler2_from_frame(frame2_from_angles);

			TS_ASSERT_EULER_EQUALS(euler2, euler2_from_frame);
			TS_ASSERT_MATRIX_EQUALS(frame2, frame2_from_angles);

			Matrix rot6, rot7, rot8, rot9;

			rot6 = z_rotation_matrix_degrees( 30. );
			rot7 = z_rotation_matrix_degrees( 250. );
			rot8 = x_rotation_matrix_degrees( 178. );
			rot9 = x_rotation_matrix_degrees( 182. );

			Matrix frame3 = rot6 * rot8 * rot7;
			Matrix frame4 = rot6 * rot9 * rot7;

			Angles euler3(frame3);
			Angles euler4(frame4);

			/// When phi and psi are "flipped" as z moves from < 180 to > 180,
			/// they are move by 180 degrees
			TS_ASSERT_DELTA( 180 + euler4.phi_degrees(), euler3.phi_degrees(), 1e-6 );
			TS_ASSERT_DELTA( 180 + euler4.psi_degrees(), 360 + euler3.psi_degrees(), 1e-6 );
		}

		void test_distance_calculation()
		{
			Matrix rot1, rot2, id;
			rot1 = y_rotation_matrix_degrees( 30. );
			rot2 = y_rotation_matrix_degrees( 45. );

			Angles euler1(rot1);
			Angles euler2(rot2);
			Angles euler_id;

			TS_ASSERT_DELTA(30., conversions::degrees(Angles::angular_distance_between(euler_id, euler1)), 1e-6);
			TS_ASSERT_DELTA(30., conversions::degrees(Angles::angular_distance_between(euler1, euler_id)), 1e-6);

			TS_ASSERT_DELTA(45., conversions::degrees(Angles::angular_distance_between(euler_id, euler2)), 1e-6);
			TS_ASSERT_DELTA(45., conversions::degrees(Angles::angular_distance_between(euler2, euler_id)), 1e-6);

			TS_ASSERT_DELTA(15., conversions::degrees(Angles::angular_distance_between(euler1, euler2)), 1e-6);
			TS_ASSERT_DELTA(15., conversions::degrees(Angles::angular_distance_between(euler2, euler1)), 1e-6);
		}
};
