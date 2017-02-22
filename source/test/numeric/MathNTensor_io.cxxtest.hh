// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathNTensor_io.cxxtest.hh
/// @brief  Test suite for numeric::MathNTensor_io
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <basic/database/open.hh>
#include <utility/tools/make_vector1.hh>
#include <numeric/MathNTensor_io.hh>
#include <numeric/types.hh>

// --------------- Test Class --------------- //

class MathNTensorIOTests : public CxxTest::TestSuite {

	public:

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //
	void test_readTensorFromFile() {

		// drawn from tests on "loop_close" 6D potential for
		// a loop of length 1; see  apps/pilot/rhiju/read_tensor.cc

		using utility::tools::make_vector1;

		std::cout << "scoring/loop_close/6D_potentials/rna/loop_01/potential.txt.gz" << std::endl;
		std::string const filename = basic::database::full_name( "scoring/loop_close/6D_potentials/rna/loop_01/potential.txt.gz" );


		numeric::MathNTensor< double, 6 > T;

		read_tensor_from_file( filename,  T );

		TS_ASSERT_EQUALS( T.n_bins().size(), 6 );
		TS_ASSERT_EQUALS( T.n_bins( 1 ), 21 );
		TS_ASSERT_EQUALS( T.n_bins( 2 ), 21 );
		TS_ASSERT_EQUALS( T.n_bins( 3 ), 21 );
		TS_ASSERT_EQUALS( T.n_bins( 4 ), 13 );
		TS_ASSERT_EQUALS( T.n_bins( 5 ), 13 );
		TS_ASSERT_EQUALS( T.n_bins( 6 ), 13 );

		utility::vector1< numeric::Size > checkbins( make_vector1( 8, 11, 11, 6, 6, 6 )  );
		TS_ASSERT_DELTA( T( make_vector1( checkbins[1]-1, checkbins[2]-1, checkbins[3]-1,
																			checkbins[4]-1, checkbins[5]-1, checkbins[6]-1 ) ),
										 2.497, 1.0e-3 );
		TS_ASSERT_DELTA( T(checkbins[1]-1, checkbins[2]-1, checkbins[3]-1,
											 checkbins[4]-1, checkbins[5]-1, checkbins[6]-1  ),
										 2.497, 1.0e-3 );

	}



};
