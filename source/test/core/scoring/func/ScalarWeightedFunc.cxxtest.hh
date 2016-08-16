// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/scoring/func/SplineFunc.cxxtest.hh
/// @brief  test suite for ScalarWeightedFunc constraints function
/// @author Stephanie DeLuca (stephanie.h.deluca@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/ScalarWeightedFunc.fwd.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

class ScalarWeightedFuncTests : public CxxTest::TestSuite{

public:
	ScalarWeightedFuncTests() {}

	// Shared initialization goes here.
	void setUp()
	{
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	void test_ScalarWeightedFunc()
	{
		static basic::Tracer TR("core.scoring.constraints.ScalarWeightedFunc.cxxtest");
		using namespace core::scoring::func;

		core::Real weight=0.0;
		FuncOP myfunc;

		ScalarWeightedFuncOP func( new ScalarWeightedFunc(weight,myfunc) );

		// Get the scoring/constraints/epr_distance_potential.histogram from mini database
		std::string epr_dist_histogram( basic::database::full_name("scoring/constraints/epr_distance_potential.histogram"));

		// Need to give some test input
		std::stringstream test_input;
		std::string line;
		test_input << "1.0\tSPLINE\tEPR_DISTANCE\t" << epr_dist_histogram << "\t0.0\t1.0\t0.5\n"
			"2.0\tBOUNDED\t-7.0\t17.0\t1.0\ttest_sd1.0\n"
			"0.5\tHARMONIC\t1.0\t1.0\n"
			"5.0\tCONSTANTFUNC\t2.0";

		// Call the read_data() function from BoundConstraint.cc to read in the above test data
		while ( getline(test_input,line) )
				{
			// use ScalarWeightedFunc::read_data(in) to read in the test input
			func->read_data(test_input);

			// Print out the stuff that ScalarWeightedFunc function reads in and calculates
			TR << "SCALARWEIGHTEDFUNC func_to_weight:  ";
			func->show_definition(TR);

		}

		float const TOLERATED_ERROR( 0.001);

		//Test if returns right value for spline
		std::stringstream in_spline;
		in_spline << "1.0\tSPLINE\tEPR_DISTANCE\t" << epr_dist_histogram << "\t0.0\t1.0\t0.5";
		func->read_data(in_spline);

		// Array containing specific x-values to test
		const size_t func_spline_sz( 33);
		core::Real func_spline[ func_spline_sz] =
			{
			-16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0,
			1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
			};

		// Array containing the correct y-values corresponding to the values in func_x[]
		const size_t func_spline_val_sz( 33);
		core::Real func_spline_val[ func_spline_val_sz] =
			{
			0.00000, 0.00000, -2.05678e-05, -0.000288659, -0.00639566, -0.0776256, -0.157987, -0.241071,
			-0.325388, -0.410867, -0.497126, -0.582517, -0.668114, -0.752582, -0.83541, -0.916787,
			-0.990608, -0.99935, -0.998547, -0.992551, -0.981298, -0.962869, -0.935373, -0.895849,
			-0.840766, -0.76187, -0.641786, -0.445389, -0.0604065, -0.00292787, -0.00020862, 0.00000, 0.00000
			};

		// Test SplineFunc for cases whose x-values are listed in func_x[] and y-values are listed in func_val[] and func_10x[]
		for ( core::Real *pos( func_spline), *pos_end( func_spline + func_spline_sz); pos != pos_end; ++pos ) {
			TS_ASSERT_DELTA( func->func( ( *pos)*-1 ), func_spline_val[ pos-func_spline], TOLERATED_ERROR);
			func->show_definition(TR);
			TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func( ( *pos)*-1 ) << "\texpect:\t" << func_spline_val[ pos-func_spline] << std::endl;
		}

		//Test Bounded
		std::stringstream in_bound;
		in_bound << "2.0\tBOUNDED\t-7.0\t17.0\t1.0\ttest_sd1.0";
		func->read_data(in_bound);

		// Array containing specific x-values to test
		const size_t func_bound_w2_sz(41);
		core::Real func_bound_w2[ func_bound_w2_sz] =
			{
			-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
			};

		// Array containing the correct y-values corresponding to the above array of x-values.  Test for sd=1
		const size_t func_bound_w2_val_sz(41);
		core::Real func_bound_w2_val[func_bound_w2_val_sz] =
			{
			128,98,72,50,32,18,8,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5
			};

		// Test BoundFunc for cases whose x-values are listed in func_x[] and y-values are listed in func_val[] and func_10x[]
		for ( core::Real *pos( func_bound_w2), *pos_end( func_bound_w2 + func_bound_w2_sz); pos != pos_end; ++pos ) {
			TS_ASSERT_DELTA( func->func(*pos), func_bound_w2_val[ pos-func_bound_w2], TOLERATED_ERROR);
			func->show_definition(TR);
			TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func(*pos) << "\texpect:\t" << func_bound_w2_val[ pos-func_bound_w2] << std::endl;
		}

		//Test HarmonicFunc
		std::stringstream in_harmonic;
		in_harmonic << "0.5\tHARMONIC\t1.0\t1.0";
		func->read_data(in_harmonic);

		// Array containing specific x-values to test
		const size_t func_harmonic_w05_sz(11);
		core::Real func_harmonic_w05[ func_harmonic_w05_sz] =
			{
			-5,-4,-3,-2,-1,0,1,2,3,4,5
			};

		// Array containing the correct y-values corresponding to the above array of x-values.  Test for sd=1
		const size_t func_harmonic_w05_val_sz(11);
		core::Real func_harmonic_w05_val[func_harmonic_w05_val_sz] =
			{
			18,12.5,8,4.5,2,0.5,0,0.5,2,4.5,8
			};

		// Test HarmonicFunc for cases whose x-values are listed in func_x[] and y-values are listed in func_val[] and func_10x[]
		for ( core::Real *pos( func_harmonic_w05), *pos_end( func_harmonic_w05 + func_harmonic_w05_sz); pos != pos_end; ++pos ) {
			TS_ASSERT_DELTA( func->func(*pos), func_harmonic_w05_val[ pos-func_harmonic_w05], TOLERATED_ERROR);
			func->show_definition(TR);
			TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func(*pos) << "\texpect:\t" << func_harmonic_w05_val[ pos-func_harmonic_w05] << std::endl;
		}

		//Test ConstantFunc
		std::stringstream in_constant;
		in_constant << "5.0\tCONSTANTFUNC\t2.0";
		func->read_data(in_constant);

		// Array containing specific x-values to test
		const size_t func_constant_w5_sz(3);
		core::Real func_constant_w5[ func_constant_w5_sz] =
			{
			-1,0,1
			};

		// Array containing the correct y-values corresponding to the above array of x-values.  Test for sd=1
		const size_t func_constant_w5_val_sz(3);
		core::Real func_constant_w5_val[func_constant_w5_val_sz] =
			{
			10,10,10
			};

		// Test HarmonicFunc for cases whose x-values are listed in func_x[] and y-values are listed in func_val[] and func_10x[]
		for ( core::Real *pos( func_constant_w5), *pos_end( func_constant_w5 + func_constant_w5_sz); pos != pos_end; ++pos ) {
			TS_ASSERT_DELTA( func->func(*pos), func_constant_w5_val[ pos-func_constant_w5], TOLERATED_ERROR);
			func->show_definition(TR);
			TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func(*pos) << "\texpect:\t" << func_constant_w5_val[ pos-func_constant_w5] << std::endl;
		}

	}// test_ScalarWeightedFunc()


	void test_serialize_ScalarWeightedFunc() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::func;

		core::Real weight=0.0;
		FuncOP myfunc;

		ScalarWeightedFuncOP func( new ScalarWeightedFunc(weight,myfunc) );

		// Get the scoring/constraints/epr_distance_potential.histogram from mini database
		std::string epr_dist_histogram( basic::database::full_name("scoring/constraints/epr_distance_potential.histogram"));

		// Need to give some test input
		std::stringstream test_input;
		std::string line;
		test_input << "1.0\tSPLINE\tEPR_DISTANCE\t" << epr_dist_histogram << "\t0.0\t1.0\t0.5\n"
			"2.0\tBOUNDED\t-7.0\t17.0\t1.0\ttest_sd1.0\n"
			"0.5\tHARMONIC\t1.0\t1.0\n"
			"5.0\tCONSTANTFUNC\t2.0";

		// Call the read_data() function from BoundConstraint.cc to read in the above test data
		while ( getline(test_input,line) ) {
			// use ScalarWeightedFunc::read_data(in) to read in the test input
			func->read_data(test_input);
		}

		FuncOP instance( func ); // serialize this through a pointer to the base class

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		FuncOP instance2; // deserialize also through a pointer to the base class
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		// make sure the deserialized base class pointer points to a ScalarWeightedFunc
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< ScalarWeightedFunc > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
