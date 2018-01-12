// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  test/numeric/MathMatrix_serialization.cxxtest.hh
/// @brief  Unit test for checking serialization functions in MathMatrx
/// @author Rebecca Alford (rfalford12@gmail.com)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <numeric/MathMatrix.hh>
#include <numeric/MathMatrix.srlz.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR("MathMatrix_serialization");

class MathMatrix_serialization : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	void test_serialization() {

		using namespace numeric;

		TS_ASSERT( true );

#ifdef SERIALIZATION


		// Setup a dummy 3x3 matrix
		MathMatrix< core::Real > dummy( 3, 3 );
		dummy(0,0) = 4; dummy(0,1) = 6; dummy(0,2) = 8;
		dummy(1,0) = 5; dummy(1,1) = 7; dummy(1,2) = 9;
		dummy(2,0) = 1; dummy(2,1) = 2; dummy(2,2) = 3;

		// Check save and load functions
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( dummy );
		}

		MathMatrix< core::Real > dummy2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( dummy2 );
		}

		// Check that both matrices have the same size
		TS_ASSERT_EQUALS( dummy.get_number_rows(), dummy2.get_number_rows() );
		TS_ASSERT_EQUALS( dummy.get_number_cols(), dummy2.get_number_cols() );

		// Check that dummy == dummy2
		for ( core::Size ii = 0; ii < dummy.get_number_rows(); ++ii ) {
			for ( core::Size jj = 0; jj < dummy.get_number_cols(); ++jj ) {
				TS_ASSERT_EQUALS( dummy(ii,jj), dummy2(ii,jj) );
			}
		}

#endif
	}

};
