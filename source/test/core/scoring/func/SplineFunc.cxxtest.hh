// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/scoring/func/SplineFunc.cxxtest.hh
/// @brief  test suite for SplineFunc constraints function
/// @author Stephanie Hirst (stephanie.j.hirst@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/scoring/func/SplineFunc.hh>
#include <core/scoring/func/SplineFunc.fwd.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>



class SplineFuncTests : public CxxTest::TestSuite {

public:
	SplineFuncTests() {};

	// Shared initialization goes here.
	void setUp()
	{
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	void test_spline_func()
	{
		static basic::Tracer TR("core.scoring.constraints.SplineFunc.cxxtest");

		using namespace core::scoring::func;

		// func points to the SplineFunc
		SplineFuncOP func( new SplineFunc());

		// Get the scoring/constraints/epr_distance_potential.histogram from mini database
		std::string epr_dist_histogram( basic::database::full_name("scoring/constraints/epr_distance_potential.histogram"));

		std::stringstream test_input;
		test_input << "EPR_DISTANCE " << epr_dist_histogram << " 0.0 1.0 0.5";

		// Call the read_data function of the SplineFunc to get the input histogram
		// Need this to read in the input KB potential
		func->read_data( test_input);

		// Print out read-in data and bounds of histogram
		TR << "exp_val:\t" << func->get_exp_val() << std::endl;
		TR << "weight:\t" << func->get_weight() << std::endl;
		TR << "bin_size:\t" << func->get_bin_size() << std::endl;
		TR << "lbx:\t" << func->get_lower_bound_x() << std::endl;
		TR << "ubx:\t" << func->get_upper_bound_x() << std::endl;
		TR << "lby:\t" << func->get_lower_bound_y() << std::endl;
		TR << "uby:\t" << func->get_upper_bound_y() << std::endl;
		TR << "lbdy:\t" << func->get_lower_bound_dy() << std::endl;
		TR << "updy:\t" << func->get_upper_bound_dy() << std::endl;

		// Array containing specific x-values to test
		const size_t func_x_sz( 33);
		core::Real func_x[ func_x_sz] =
		{
			-16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0,
			1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
		};

		// Array containing the correct y-values corresponding to the values in func_x[]
		const size_t func_val_sz( 33);
		core::Real func_val[ func_val_sz] =
		{
			0.00000, 0.00000, -2.05678e-05, -0.000288659, -0.00639566, -0.0776256, -0.157987, -0.241071,
			-0.325388, -0.410867, -0.497126, -0.582517, -0.668114, -0.752582, -0.83541, -0.916787,
			-0.990608, -0.99935, -0.998547, -0.992551, -0.981298, -0.962869, -0.935373, -0.895849,
			-0.840766, -0.76187, -0.641786, -0.445389, -0.0604065, -0.00292787, -0.00020862, 0.00000, 0.00000
		};

		// Array containing the correct y-values * 10 corresponding to the values in func_x[]
		const size_t func_10_sz( 33);
		if ( func_x_sz == func_val_sz && func_val_sz == func_10_sz)
		{
			core::Real func_10[ func_10_sz];
			for( core::Real *pos( func_val), *pos_end( func_val + func_val_sz); pos != pos_end; ++pos)
			{
				// Multiply each value in func_val[] by 10 and put into func_10x
				func_10[ pos-func_val] = ( *pos)*10;
			}

			// Test SplineFunc for cases whose x-values are listed in func_x[] and y-values are listed in func_val[] and func_10x[]
			for( core::Real *pos( func_x), *pos_end( func_x + func_x_sz); pos != pos_end; ++pos)
			{
				float const TOLERATED_ERROR( 0.001);
				TS_ASSERT_DELTA( func->func( ( *pos)*-1 ), func_val[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func( ( *pos)*-1 ) << "\texpect:\t" << func_val[ pos-func_x] << std::endl;
				TS_ASSERT_DELTA( ( func->func( ( *pos)*-1 ))*10, func_10[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\t10func(r) =\t" << ( func->func( ( *pos)*-1 ))*10 << "\texpect:\t" << func_10[ pos-func_x] << std::endl;
			}

//			// Print out a curve with a resolution of 0.1 to see if matches with input KB potential
//			for(core::Real r( -20.0); r <= 20.0; r += 0.1)
//			{
//				// Because SplineFunc::func calculates potential from diff = exp_val_ - r, need to reverse the sign
//				TR << "r =\t" << r << "\tfunc(r) =\t" << func->func( r*-1 ) << "\t10func(r) =\t" << ( func->func( r* -1 ))*10 << std::endl;
//			}

		} // if arrays' size are equal
	} // test_spline_func()
}; // SplineFunc test suite
