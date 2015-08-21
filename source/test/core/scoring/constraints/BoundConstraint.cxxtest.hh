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

/// @file   core/scoring/constraints/BoundConstraint.cxxtest.hh
/// @brief  test suite for BoundFunc constraints function
/// @author Stephanie H. DeLuca (stephanie.h.deluca@vanderbilt.edu)

//Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.fwd.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <sstream>

//Auto Headers
//#include <utility/vector1.hh>

class BoundConstraintTests : public CxxTest::TestSuite {

public:
	BoundConstraintTests() {}

	// Shared initialization goes here.
	void setUp()
	{
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	void test_bound_func()
	{

		static basic::Tracer TR("core.scoring.constraints.BoundConstraint.cxxtest");
		using namespace core::scoring::constraints;

		core::Real lb=0;
		core::Real ub=0;
		core::Real sd=1.0;
		std::string type="";
		core::Real rswitch=0.5;

		// func points to the BoundConstraint (Func)
		BoundFuncOP func( new BoundFunc(lb,ub,sd,rswitch,type) );

		// Need to give some test input
		std::stringstream test_input;
		std::string line;
		test_input << "-7.0\t17.0\t1.0\ttest_sd1.0\n"
			"-7.0\t17.0\t0.5\ttest_sd0.5\n"
			"-7.0\t17.0\t2.5\ttest_sd2.5\n"
			"-7.0\t17.0\t5.0\ttest_sd5.0";

		// Call the read_data() function from BoundConstraint.cc to read in the above test data
		while ( getline(test_input,line) )
				{
			func->read_data(test_input);

			// Print out the stuff that BoundConstraint function reads in and calculates
			TR << "lower_bound:\t" << func->lb() << std::endl;
			TR << "upper_bound:\t" << func->ub() << std::endl;
			TR << "sd:\t" << func->sd() << std::endl;
			TR << "rswitch:\t" << func->rswitch() << std::endl;
			TR << "type:\t" << func->type() << std::endl;

			if ( rswitch!=0.5 || func->rswitch()!=0.5 ) {
				TS_ASSERT(rswitch==0.5);
				TS_ASSERT(func->rswitch()==0.5);
			}

		}

		// Array containing specific x-values to test
		const size_t func_x_sz(41);
		core::Real func_x[ func_x_sz] =
			{
			-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
			};

		// Array containing the correct y-values corresponding to the above array of x-values.  Test for sd=1
		const size_t func_sd1_w1_sz(41);
		core::Real func_sd1_w1[func_sd1_w1_sz] =
			{
			64,49,36,25,16,9,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.75,1.75,2.75,3.75,4.75,5.75,6.75,7.75
			};

		// Now test to for sd=0.5, 2.5, and 5.0
		const size_t func_sd0point5_w1_sz(41);
		core::Real func_sd0point5_w1[func_sd0point5_w1_sz]=
			{
			256,196,144,100,64,36,16,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.75,3.75,5.75,7.75,9.75,11.75,13.75,15.75
			};

		const size_t func_sd2point5_w1_sz(41);
		core::Real func_sd2point5_w1[func_sd2point5_w1_sz]=
			{
			10.24,7.84,5.76,4,2.56,1.44,0.64,0.16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.16,0.55,0.95,1.35,1.75,2.15,2.55,2.95
			};

		const size_t func_sd5_w1_sz(41);
		core::Real func_sd5_w1[func_sd5_w1_sz]=

			{
			2.56,1.96,1.44,1,0.64,0.36,0.16,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.16,0.35,0.55,0.75,0.95,1.15,1.35
			};

		// Array containing the correct y-values if weight=5.0, 10.0, 15.0, 20.0
		const size_t func_5_sz(41);
		core::Real func_5[func_5_sz];

		const size_t func_10_sz(41);
		core::Real func_10[func_10_sz];

		const size_t func_15_sz(41);
		core::Real func_15[func_15_sz];

		const size_t func_20_sz(41);
		core::Real func_20[func_20_sz];

		// If everything is the right size, continue
		if ( func_x_sz == func_sd1_w1_sz && func_sd1_w1_sz == func_10_sz && func_sd1_w1_sz == func_15_sz &&
				func_sd1_w1_sz == func_20_sz && func_sd0point5_w1_sz==func_sd1_w1_sz && func_sd2point5_w1_sz==func_sd1_w1_sz &&
				func_sd5_w1_sz==func_sd1_w1_sz ) {
			// fill the func vectors with expected values
			for ( core::Real *pos(func_sd1_w1),*pos_end(func_sd1_w1+func_sd1_w1_sz); pos!=pos_end; ++pos ) {
				func_5[pos-func_sd1_w1] = (*pos)*5;
				func_10[pos-func_sd1_w1] = (*pos)*10;
				func_15[pos-func_sd1_w1] = (*pos)*15;
				func_20[pos-func_sd1_w1] = (*pos)*20;
			}

			// fill func_sd0point5_w1,func_sd2point5_w1, and func_sd5_w1
			for ( core::Real *pos(func_sd0point5_w1),*pos_end(func_sd0point5_w1+func_sd0point5_w1_sz);
					pos!=pos_end; ++pos ) {
				func_sd0point5_w1[pos-func_sd0point5_w1] = (*pos);
			}

			for ( core::Real *pos(func_sd2point5_w1),*pos_end(func_sd2point5_w1+func_sd2point5_w1_sz); pos!=pos_end; ++pos ) {
				func_sd2point5_w1[pos-func_sd2point5_w1] = (*pos);
			}

			for ( core::Real *pos(func_sd5_w1),*pos_end(func_sd5_w1+func_sd5_w1_sz); pos!=pos_end; ++pos ) {
				func_sd5_w1[pos-func_sd5_w1] = (*pos);
			}


			// Test BoundFunc for cases whose x-values are listed in func_x[] and y-values are listed in func_val[] and func_10x[]
			for ( core::Real *pos( func_x), *pos_end( func_x + func_x_sz); pos != pos_end; ++pos ) {
				float const TOLERATED_ERROR( 0.001);
				std::stringstream in_sd1;
				in_sd1 << "-7.0\t17.0\t1.0\ttest_sd1.0";
				func->read_data(in_sd1);
				TS_ASSERT_DELTA(func->sd(),1.0,TOLERATED_ERROR);

				TR << "testing sd=1,w=1" << std::endl;
				TS_ASSERT_DELTA( func->func(*pos), func_sd1_w1[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func(*pos) << "\texpect:\t" << func_sd1_w1[ pos-func_x] << std::endl;

				TR << "testing sd=1,w=5" << std::endl;
				TS_ASSERT_DELTA(func->func(*pos)*5,func_5[pos-func_x],TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\t5func(r) =\t" << ( func->func(*pos))*5 << "\texpect:\t" << func_5[ pos-func_x] << std::endl;

				TR << "testing sd=1,w=10" << std::endl;
				TS_ASSERT_DELTA( ( func->func( ( *pos)))*10, func_10[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\t10func(r) =\t" << ( func->func(*pos))*10 << "\texpect:\t" << func_10[ pos-func_x] << std::endl;

				TR << "testing sd=1,w=15" << std::endl;
				TS_ASSERT_DELTA( ( func->func( ( *pos)))*15, func_15[pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\t15func(r) =\t" << ( func->func(*pos))*15 << "\texpect:\t" << func_15[ pos-func_x] << std::endl;

				TR << "testing sd=1,w=20" << std::endl;
				TS_ASSERT_DELTA( ( func->func( ( *pos)))*20, func_20[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\t20func(r) =\t" << ( func->func(*pos))*20 << "\texpect:\t" << func_20[ pos-func_x] << std::endl;

				std::stringstream in_sd0point5;
				in_sd0point5 << "-7.0\t17.0\t0.5\ttest_sd0.5";
				func->read_data(in_sd0point5);
				TS_ASSERT_DELTA(func->sd(),0.5,TOLERATED_ERROR);
				TR << "testing sd=0.5,w=1" << std::endl;
				TS_ASSERT_DELTA( func->func(*pos), func_sd0point5_w1[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func(*pos) << "\texpect:\t" << func_sd0point5_w1[ pos-func_x] << std::endl;

				std::stringstream in_sd2point5;
				in_sd2point5 << "-7.0\t17.0\t2.5\ttest_sd2.5";
				func->read_data(in_sd2point5);
				TS_ASSERT_DELTA(func->sd(),2.5,TOLERATED_ERROR);
				TR << "testing sd=2.5,w=1" << std::endl;
				TS_ASSERT_DELTA( func->func(*pos), func_sd2point5_w1[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func(*pos) << "\texpect:\t" << func_sd2point5_w1[ pos-func_x] << std::endl;

				std::stringstream in_sd5;
				in_sd5 << "-7.0\t17.0\t5.0\ttest_sd5";
				func->read_data(in_sd5);
				TS_ASSERT_DELTA(func->sd(),5.0,TOLERATED_ERROR);
				TR << "testing sd=5,w=1" << std::endl;
				TS_ASSERT_DELTA( func->func(*pos), func_sd5_w1[ pos-func_x], TOLERATED_ERROR);
				TR << "r =\t" << *pos << "\tfunc(r) =\t" << func->func(*pos) << "\texpect:\t" << func_sd5_w1[ pos-func_x] << std::endl;

			} // for( core::Real *pos( func_x), *pos_end( func_x + func_x_sz); pos != pos_end; ++pos)

			//    // Print out a curve with a resolution of 0.1 to see if matches with what we expect
			//    for(core::Real r( -25.0); r <= 25.0; r += 0.1)
			//    {
			//     TR << "sd = " << func->sd() << "r = " << r << "\tfunc(r) = " << func->func(r) << std::endl;
			//    }

		} // if(func_x_sz == func_sd1_w1_sz)

	}// test_bound_func()

	//at SVN 49611, this function was broken because it used setprecision() without ever *resetting* the precision.  Oops.
	void test_show_definition(){

		//Make a quick BoundConstraint
		core::Real lb=0;
		core::Real ub=0;
		core::Real sd=1.0;
		std::string type="";
		core::Real rswitch=0.5;
		core::scoring::constraints::BoundFuncOP func( new core::scoring::constraints::BoundFunc(lb,ub,sd,rswitch,type) );

		std::stringstream test_stream, compare_stream;

		//these are random (well, convenience) numbers to test the precision problem in BoundFunc
		//note 4.9938572 is shortened to 4.99386 because the default precision is apparently 5
		test_stream << 4.9938572 << " 5.6732247" << std::endl;
		func->show_definition(test_stream);
		test_stream << 4.9938572 << " 5.6732247" << std::endl;
		compare_stream << test_stream.str() << std::endl;
		//std::cout << compare_stream.str();
		//std::cout << "blank" << std::endl;
		TS_ASSERT(compare_stream.str() == "4.99386 5.6732247\nBOUNDED       0       0   1 \n4.99386 5.6732247\n\n");

	}
}; // BoundedFunc test suite
