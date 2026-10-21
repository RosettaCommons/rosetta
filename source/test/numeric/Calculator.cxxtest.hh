// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/CalculatorTest.cxxtest.hh
/// @brief
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <numeric/Calculator.hh>

#include <map>

using namespace numeric;

class CalculatorTest : public CxxTest::TestSuite {
public:

	numeric::Real delta_;

	CalculatorTest() {};

	// Shared initialization goes here.
	void setUp() {
		delta_ = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {}

	void test_simple() {
		numeric::Calculator calc("3.14    ");

		std::map< std::string, numeric::Real > vars;
		numeric::Real out = 9999;

		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 3.14, delta_ );
	}

	void test_expression() {
		numeric::Calculator calc(" 1 + 2 * 3 - 6 / ( 3 + 9 ) + -(2) ");

		std::map< std::string, numeric::Real > vars;
		numeric::Real out = 9999;

		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 4.5 , delta_ );
	}

	void test_exponent() {
		numeric::Calculator calc("1.5 * 2 ^ 3 ^ 2 / 2"); // 1.5 * ( 2^(3^2) ) / 2

		std::map< std::string, numeric::Real > vars;
		numeric::Real out = 9999;

		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 384 , delta_ );
	}

	void test_variables() {
		numeric::Calculator calc(" alpha + beta ");

		std::map< std::string, numeric::Real > vars;
		vars["alpha"] = 4;
		vars["beta"] = 5;
		numeric::Real out = 9999;

		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 9 , delta_ );
	}

	void test_assign() {
		numeric::Calculator calc(" alpha = 1 + 3; beta = 3/4; alpha * beta ");

		std::map< std::string, numeric::Real > vars;
		numeric::Real out = 9999;

		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 3 , delta_ );
	}

	void test_functions() {
		std::map< std::string, numeric::Real > vars;
		vars["cosy"] = 1.5;
		vars["logs"] = 8;
		numeric::Real out = 9999;

		numeric::Calculator calc("log(100)");

		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 2 , delta_ );

		numeric::Calculator calc2(" cos(SIN( d2r(-180)) ) + ln(exp(log(10)+2)) + sqrt(9)*LOG2(logs) + cosy*(2-1) + log(27,3) + mean(1, 1.5, 3, 4.5) + ABS(min(-0.2, 3, 0.6))"); // 1 + 3 + 3*3 + 1.5 + 3 + 2.5 + 0.2

		out = 9999;
		TS_ASSERT( ! calc2.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 20.2 , delta_ );
	}

	void test_sinpass() {

		numeric::Calculator calc("sinpass(x, 30.0, 50.0, 10.0, 16.0)");
		std::map< std::string, numeric::Real > vars;
		numeric::Real out = 9999;

		//for ( core::Size i(0); i <= 1000; i++ ) {
		// vars["x"] = i * 0.1;
		// calc.compute( vars, out );
		// std::cout << vars["x"] << "   " << out << std::endl;
		//}

		vars["x"] = 30;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 1.0 , delta_ );

		vars["x"] = 43.2;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 1.0 , delta_ );

		vars["x"] = 50;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 1.0 , delta_ );

		vars["x"] = 30-0.00001;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 1.0 , delta_ ); // Not exactly, but within delta

		vars["x"] = 50+0.00001;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 1.0 , delta_ ); // Not exactly, but within delta

		vars["x"] = 20;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.0 , delta_ );

		vars["x"] = 66;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.0 , delta_ );

		vars["x"] = 20+0.00001;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.0 , delta_ ); // Not exactly, but within delta

		vars["x"] = 66-0.00001;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.0 , delta_ ); // Not exactly, but within delta

		vars["x"] = -135.7;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.0 , delta_ );

		vars["x"] = 2673.4;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.0 , delta_ );

		vars["x"] = 25;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.5 , delta_ );

		vars["x"] = 58;
		TS_ASSERT( ! calc.compute( vars, out ) );
		TS_ASSERT_DELTA( out, 0.5 , delta_ );

	}
};
