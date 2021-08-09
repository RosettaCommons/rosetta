// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  utility/temp/BinaryUtilTests.cxxtest.hh
/// @brief  Unit tests for the encode6bit and decodet6bit functions.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <utility/Binary_Util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("BinaryUtilTests");


class BinaryUtilTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}

	/// @brief Ensure that we can write to a chunk of memory without buffer overruns.
	void test_decode6bit_avoid_buffer_overrun() {
		// The following deliberately uses doubles and not Reals.
		double myvec[20];
		for ( core::Size i(0); i<20; ++i ) {
			myvec[i] = 6.34;
		}

		std::string compressed;


		//Encode the data above.
		std::string mystring;
		utility::encode6bit( reinterpret_cast< unsigned char const * >( &myvec[0] ), 20*sizeof(double), mystring );
		TR << "String after encoding: " << mystring << std::endl;

		//Append some characters.
		mystring += "TrailingCharacters";
		TR << "String after appending nonsense: " << mystring << std::endl;

		//Create a new vector to receive the data that are now longer.
		double myvec2[50];
		for ( core::Size i(0); i<50; ++i ) {
			myvec2[i] = 7.11;
		}

		//Decode the data into the new vector, stopping when we run out of data.
		utility::decode6bit( reinterpret_cast< unsigned char * >( &myvec2[0] ), mystring, 20*sizeof(double) );

		//Check that we've stopped at the end of the data and not overrun.
		TR << "Decoded vector: ";
		for ( core::Size i(0); i<50; ++i ) {
			TR << myvec2[i];
			if ( i<49 ) {
				TR << ", ";
			} else {
				TR << std::endl;
			}
			TS_ASSERT_DELTA( myvec2[i], (i < 20 ? 6.34 : 7.11), 0.001 );
		}
	}

};
