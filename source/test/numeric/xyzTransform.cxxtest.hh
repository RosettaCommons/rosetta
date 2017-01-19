// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyzTransform.cxxtest.hh
/// @brief  test suite for numeric::xyzTransform
/// @author Darwin Y. Fu


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>
#include <numeric/xyzTransform.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/numbers.hh>
#include <iostream>
#include <core/types.hh>


// --------------- Test Class --------------- //

class XYZTransformTests : public CxxTest::TestSuite {

	public:

	utility::vector1<core::PointPosition> coords_;
	core::PointPosition center_;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		//Triangular Pyramid centered at (1.25, 1.25, 1.25)
		coords_.push_back(core::PointPosition(2,1,1));
		coords_.push_back(core::PointPosition(1,2,1));
		coords_.push_back(core::PointPosition(1,1,2));
		coords_.push_back(core::PointPosition(1,1,1));
		center_ = core::PointPosition(1.25,1.25,1.25);

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief 3D Transformation of coordinates with translation and rotation about centroid
	void test_3D_transformation() {

		//Setup
		numeric::xyzMatrix<core::Real> rotation(
				numeric::z_rotation_matrix_degrees( 115.0 ) * (
					numeric::y_rotation_matrix_degrees( 85.0 ) *
					numeric::x_rotation_matrix_degrees( 5.0 ) ));

		core::Vector translation(1.5,2.5,3.0);

		//Correct Values from manual calculation:
		utility::vector1<core::PointPosition> correct_values;
		correct_values.push_back(core::PointPosition(3.04237,3.66076,3.47925));
		correct_values.push_back(core::PointPosition(2.13965,3.23945,4.48304));
		correct_values.push_back(core::PointPosition(2.73878,4.51803,4.56227));
		correct_values.push_back(core::PointPosition(3.0792,3.58177,4.47544));

		utility::vector1<core::PointPosition> test_coords = coords_;

		core::PointPosition new_center = center_ + translation;

		//Get Transformer and perform transformation
		numeric::xyzTransform<core::Real> transformer(numeric::xyzTransform<core::Real>::rot(rotation,center_,new_center));

		for(utility::vector1<core::PointPosition>::iterator it = test_coords.begin(); it != test_coords.end(); ++it)
		{
			*it = transformer*(*it);
		}

		for ( core::Size atom_index = 1; atom_index <= test_coords.size(); ++atom_index ) {
			//std::cout << atom_index << "," << test_coords[atom_index].x() << "," << test_coords[atom_index].y() << "," << test_coords[atom_index].z() << std::endl;
			TS_ASSERT_DELTA(correct_values[atom_index],test_coords[atom_index],0.001);
		}

	}

};


