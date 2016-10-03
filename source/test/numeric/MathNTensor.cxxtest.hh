// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathNTensor.cxxtest.hh
/// @brief  Test suite for numeric::MathNTensor
/// @author Vikram K. Mulligan (vmullig@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/MathNTensor.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/types.hh>

// --------------- Test Class --------------- //

class MathNTensorTests : public CxxTest::TestSuite {

	public:

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	/// @brief Test conversion to a MathMatrix.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_convert_to_MathMatrix() {
		using namespace numeric;
		
		utility::vector1< Size > dimensions(2);
		dimensions[1] = 4; dimensions[2] = 2;
		
		//Make a MathNTensor and fill it with 1..8.
		MathNTensor< Size, 2 > mytensor( dimensions );
		Size count=1;
		for(Size j=0; j<2; ++j) {
			for(Size i=0; i<4; ++i) {
				utility::vector1 <Size> coords(2);
				coords[1]=i; coords[2]=j;
				mytensor( coords ) = count;
				++count;
			}
		}
		
		//Convert to a MathMatrix and verify contents:
		MathMatrix <Size> mymatrix( mytensor.get_mathmatrix() );
		Size count2=1;
		for(Size j=0; j<2; ++j) {
			for(Size i=0; i<4; ++i) {
				TS_ASSERT_EQUALS( mymatrix(i,j), count2 );
				++count2;
			}
		}
	}



};
