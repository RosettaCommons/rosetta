// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/interpolation/spline
/// @brief  test suite for numeric::interpolation::spline::Bicubic_spline
/// @author Steven Combs (steven.combs@vanderbilt.edu)
/// This tests the functions that are in the bicubic spline class except for
/// the e_periodic steps.


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/constants.hh>
#include <numeric/MathTensor.hh>
#include <numeric/MathNTensor.hh>
#include <numeric/interpolation/spline/PolycubicSpline.hh>
#include <numeric/interpolation/spline/PolycubicSpline.tmpl.hh>
#include <numeric/interpolation/spline/TricubicSpline.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>
#include <numeric/interpolation/spline/CubicSpline.fwd.hh>
#include <numeric/MathVector_operations.hh>
#include <utility/fixedsizearray1.hh>


// --------------- Test Class --------------- //


class PolycubicSplineTests : public CxxTest::TestSuite {


public:
	//shared data


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp()
	{
	}


	// Shared finalization goes here.
	void tearDown() {

	}
	
	void test_polycubic_spline_4_data_access() {

		using namespace numeric;
		using namespace numeric::interpolation::spline;
		const double values[] = {
			/// Superlayer 1
			//layer 0 (x0)
			//z->
			26,3, 1, 2, 1, 3, //y
			6, 3, 8, 2, 7, 8,
			3, 4, 2, 1, 2, 5,
			30,0, 2, 4, 6, 3,
			4, 3, 3, 4, 11,5,
			8, 5, 2, 0, 2, 2,
			//layer 1 (x1)
			//z->
			2, 0, 0, 1, 2, 0, //y
			2, 1, 1, 3, 3, 0,
			1, 0, 0, 1, 0, 0,
			4, 0, 0, 0, 1, 3,
			0, 5, 0, 0, 1, 1,
			1, 0, 0, 3, 0, 0,
			//layer 2 (x2)
			//z->
			1, 0, 0, 0, 0, 0, //y
			1, 0, 0, 1, 2, 3,
			0, 0, 0, 1, 0, 0,
			76,0, 0, 0, 0, 0,
			0, 3, 2, 2, 0, 0,
			0, 0, 0, 0, 0, 0,
			//layer 3 (x3)
			//z->
			476, 107, 19, 2, 2, 0, //y
			0,     0,  0, 0, 0, 0,
			0,    1, 9, 3, 10, 14,
			447, 14, 13, 4, 9,  5,
			1,  3,  0, 0, 0,  0,
			0, 0, 1, 5, 14, 104,
			//layer 4 (x4)
			//z->
			1, 1, 4, 2, 1, 2, //y
			2, 1, 1, 0, 1, 0,
			2, 1, 0, 0, 0, 1,
			35,0, 0, 0, 4, 0,
			0, 1, 4, 3, 0, 0,
			0, 2, 2, 1, 1, 1,

			/// Superlayer 2
			//layer 0 (x0)
			//z->
			26,3, 1, 2, 1, 3, //y
			6, 3, 8, 2, 7, 8,
			3, 4, 2, 1, 2, 5,
			30,0, 2, 4, 6, 3,
			4, 3, 3, 4, 11,5,
			8, 5, 2, 0, 2, 2,
			//layer 1 (x1)
			//z->
			2, 0, 0, 1, 2, 0, //y
			2, 1, 1, 3, 3, 0,
			1, 0, 0, 1, 0, 0,
			4, 0, 0, 0, 1, 3,
			0, 5, 0, 0, 1, 1,
			1, 0, 0, 3, 0, 0,
			//layer 2 (x2)
			//z->
			1, 0, 0, 0, 0, 0, //y
			1, 0, 0, 1, 2, 3,
			0, 0, 0, 1, 0, 0,
			76,0, 0, 0, 0, 0,
			0, 3, 2, 2, 0, 0,
			0, 0, 0, 0, 0, 0,
			//layer 3 (x3)
			//z->
			476, 107, 19, 2, 2, 0, //y
			0,     0,  0, 0, 0, 0,
			0,    1, 9, 3, 10, 14,
			447, 14, 13, 4, 9,  5,
			1,  3,  0, 0, 0,  0,
			0, 0, 1, 5, 14, 104,
			//layer 4 (x4)
			//z->
			1, 1, 4, 2, 1, 2, //y
			2, 1, 1, 0, 1, 0,
			2, 1, 0, 0, 0, 1,
			35,0, 0, 0, 4, 0,
			0, 1, 4, 3, 0, 0,
			0, 2, 2, 1, 1, 1,

			/// Superlayer 3
			//layer 0 (x0)
			//z->
			26,3, 1, 2, 1, 3, //y
			6, 3, 8, 2, 7, 8,
			3, 4, 2, 1, 2, 5,
			30,0, 2, 4, 6, 3,
			4, 3, 3, 4, 11,5,
			8, 5, 2, 0, 2, 2,
			//layer 1 (x1)
			//z->
			2, 0, 0, 1, 2, 0, //y
			2, 1, 1, 3, 3, 0,
			1, 0, 0, 1, 0, 0,
			4, 0, 0, 0, 1, 3,
			0, 5, 0, 0, 1, 1,
			1, 0, 0, 3, 0, 0,
			//layer 2 (x2)
			//z->
			1, 0, 0, 0, 0, 0, //y
			1, 0, 0, 1, 2, 3,
			0, 0, 0, 1, 0, 0,
			76,0, 0, 0, 0, 0,
			0, 3, 2, 2, 0, 0,
			0, 0, 0, 0, 0, 0,
			//layer 3 (x3)
			//z->
			476, 107, 19, 2, 2, 0, //y
			0,     0,  0, 0, 0, 0,
			0,    1, 9, 3, 10, 14,
			447, 14, 13, 4, 9,  5,
			1,  3,  0, 0, 0,  0,
			0, 0, 1, 5, 14, 104,
			//layer 4 (x4)
			//z->
			1, 1, 4, 2, 1, 2, //y
			2, 1, 1, 0, 1, 0,
			2, 1, 0, 0, 0, 1,
			35,0, 0, 0, 4, 0,
			0, 1, 4, 3, 0, 0,
			0, 2, 2, 1, 1, 1
      };
		
		//the tensor describes layers, rows, columns, source of input
		utility::fixedsizearray1< Size, 4 > bins;
		bins[1] = 3;
		bins[2] = 5;
		bins[3] = 6;
		bins[4] = 6;
		const MathNTensor< Real, 4 > input_tensor( bins, values);
		
		utility::fixedsizearray1< BorderFlag, 4 > border;
		border[1] = e_Natural;
		border[2] = e_Natural;
		border[3] = e_Natural;
		border[4] = e_Periodic;
		
		//these vectors are used to input the starting points and the
		//grid width delta of every dimension (x, y, z) into the spline
		utility::fixedsizearray1< double, 4 > const start( 10 );
		utility::fixedsizearray1< double, 4 > const delta( 10 );
		
		utility::fixedsizearray1< bool, 4 > const lin_cont( true );
		
		//this vector controls the behavior of the spline at the beginning and
		//end of every dimension, only has impact for BorderFlag FIRSTDER
		
		//every pair describes the value of the first order derivative at
		//start and end
		utility::fixedsizearray1< std::pair< double, double>, 4 > const first_be( std::pair< double, double>( 10, 10) );
		
		
		// natural means, the second order derivative at beginning and end is 0
		// this leads to a rather slow changing first order derivative meaning
		// "nearly linear" behavior of the spline
		// in z-direction the spline is made periodic
		PolycubicSpline< 4 > naturalspline;
		naturalspline.train( border, start, delta, input_tensor, lin_cont, first_be);
		
		for( double z( 20); z < 60; z += 4)
		{
			utility::fixedsizearray1< Real, 4 > vals(10);
			utility::fixedsizearray1< Real, 4 > valsp60(10);
			vals[4] = z;
			valsp60[4] = z + 60;
			
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline.F( valsp60 ), 1e-12 );
		}
	}
	
	void test_polycubic_spline_3_data_access() {

		using namespace numeric;
		using namespace numeric::interpolation::spline;
		const double values[] = {
			//layer 0 (x0)
			//z->
			26,3, 1, 2, 1, 3, //y
			6, 3, 8, 2, 7, 8,
			3, 4, 2, 1, 2, 5,
			30,0, 2, 4, 6, 3,
			4, 3, 3, 4, 11,5,
			8, 5, 2, 0, 2, 2,
			//layer 1 (x1)
			//z->
			2, 0, 0, 1, 2, 0, //y
			2, 1, 1, 3, 3, 0,
			1, 0, 0, 1, 0, 0,
			4, 0, 0, 0, 1, 3,
			0, 5, 0, 0, 1, 1,
			1, 0, 0, 3, 0, 0,
			//layer 2 (x2)
			//z->
			1, 0, 0, 0, 0, 0, //y
			1, 0, 0, 1, 2, 3,
			0, 0, 0, 1, 0, 0,
			76,0, 0, 0, 0, 0,
			0, 3, 2, 2, 0, 0,
			0, 0, 0, 0, 0, 0,
			//layer 3 (x3)
			//z->
			476, 107, 19, 2, 2, 0, //y
			0,     0,  0, 0, 0, 0,
			0,    1, 9, 3, 10, 14,
			447, 14, 13, 4, 9,  5,
			1,  3,  0, 0, 0,  0,
			0, 0, 1, 5, 14, 104,
			//layer 4 (x4)
			//z->
			1, 1, 4, 2, 1, 2, //y
			2, 1, 1, 0, 1, 0,
			2, 1, 0, 0, 0, 1,
			35,0, 0, 0, 4, 0,
			0, 1, 4, 3, 0, 0,
			0, 2, 2, 1, 1, 1
      };
		
		utility::fixedsizearray1< Size, 3 > dims;
		dims[1] = 5; dims[2] = 6; dims[3] = 6; 
		//the tensor describes layers, rows, columns, source of input
		const MathNTensor< Real, 3 > input_tensor( dims, values);
		
		utility::fixedsizearray1< BorderFlag, 3 > border( e_Natural );
		border[3] = e_Periodic;
		
		//these vectors are used to input the starting points and the
		//grid width delta of every dimension (x, y, z) into the spline
		utility::fixedsizearray1< double, 3> start( 10 );
		utility::fixedsizearray1< double, 3> delta( 10 );
		utility::fixedsizearray1< bool, 3> lin_cont( true );
		utility::fixedsizearray1< std::pair< double, double>, 3> first_be( std::pair< double, double>( 10, 10) );
		
		// natural means, the second order derivative at beginning and end is 0
		// this leads to a rather slow changing first order derivative meaning
		// "nearly linear" behavior of the spline
		// in z-direction the spline is made periodic
		PolycubicSpline< 3 > naturalspline;
		naturalspline.train( border, start, delta, input_tensor, lin_cont, first_be);
		
		//the next line generates a class from a trained spline that can give back F(x), dF(x) and FdF(x)
		//naturalspline.WriteCodeInFile("natural_trained", GetExampleSourceCodePath());
		
		//std::cout << "To show continuous behavior at the end of a cell" << std::endl;
		//std::cout << " x        F(x, 10, 10)       F_x(x, 10, 10)" << std::endl;
		for( double x( 19.9); x < 20.1; x += 0.01)
		{
			//std::cout << x << " " << naturalspline.F( x, 10, 10) << " " << naturalspline.dFdx( x, 10, 10) << std::endl;
		}
		
		//std::cout << "Behavior at the end of the defined region" << std::endl;
		for( double x( 49); x < 50.5; x += 0.1)
		{
			//std::cout << x << " ";
			//std::cout << naturalspline.F( x, 10, 10) << " ";
			//std::cout << naturalspline.dFdx( x, 10, 10) << std::endl;
		}
		//std::cout << "To show periodic behavior in z-direction" << std::endl;
		//std::cout << " z        F(10, 10, z)       F(10, 10, z+60)" << std::endl;
		for( double z( 20); z < 60; z += 4)
		{
			//std::cout <<  z <<d " ";
			//std::cout << naturalspline.F( 10, 10, z) << " ";
			//std::cout <<  naturalspline.F( 10, 10, z+60) << std::endl;
			utility::fixedsizearray1< Real, 3 > vals( 10.0 );
			utility::fixedsizearray1< Real, 3 > valsp60( 10.0 );
			vals[3] = z;
			valsp60[3] = z+60;
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline.F( valsp60 ), 1e-12 );
		}
		
		//this example describes a function f(x, y, z), where (x, y, z) resembles a helix
		//the spline is periodic in every direction, but this is not shown here
		PolycubicSpline< 3 > ts;
		
		//std::cout <<  "Training new purely periodic spline" << std::endl;
		
		border[ 1 ] = e_Periodic;
		border[ 2 ] = e_Periodic;
		border[ 3 ] = e_Periodic;

		ts.train( border, start, delta, input_tensor, lin_cont, first_be);
		
		//last function example
		//std::cout <<  "Function values along a helix" << std::endl;
		//std::cout << "x F( 10+x, 10+10*std::cos( pi / 18 * x), 10+10*std::sin( pi / 18 * x))" << std::endl;
		for( double x( 0); x <= 36; ++x) {
			//std::cout <<  x << " ";
			//std::cout << ts.F( 10 + x, 10 + 10 * std::cos( numeric::constants::d::pi / 18 * x), 10 + 10 * std::sin( numeric::constants::d::pi / 18 * x)) << std::endl;
		}
	}

	void test_tricubic_polycubic_3_identical() {
		using namespace numeric;
		using namespace numeric::interpolation;
		using namespace numeric::interpolation::spline;

		const double values[] = {
			//layer 0 (x0)
			//z->
			 26,   3,   1,   2,   1,   3, //y
			  6,   3,   8,   2,   7,   8,
			  3,   4,   2,   1,   2,   5,
			 30,   0,   2,   4,   6,   3,
			  4,   3,   3,   4,  11,   5,
			  8,   5,   2,   0,   2,   2,
			//layer 1 (x1)
			//z->
			  2,   0,   0,   1,   2,   0, //y
			  2,   1,   1,   3,   3,   0,
			  1,   0,   0,   1,   0,   0,
			  4,   0,   0,   0,   1,   3,
			  0,   5,   0,   0,   1,   1,
			  1,   0,   0,   3,   0,   0,
			//layer 2 (x2)
			//z->
			  1,   0,   0,   0,   0,   0, //y
			  1,   0,   0,   1,   2,   3,
			  0,   0,   0,   1,   0,   0,
			 76,   0,   0,   0,   0,   0,
			  0,   3,   2,   2,   0,   0,
			  0,   0,   0,   0,   0,   0,
			//layer 3 (x3)
			//z->
			476, 107,  19,   2,   2,   0, //y
			  0,   0,   0,   0,   0,   0,
			  0,   1,   9,   3,  10,  14,
			447,  14,  13,   4,   9,   5,
			  1,   3,   0,   0,   0,   0,
			  0,   0,   1,   5,  14, 104,
			//layer 4 (x4)
			//z->
			  1,   1,   4,   2,   1,   2, //y
			  2,   1,   1,   0,   1,   0,
			  2,   1,   0,   0,   0,   1,
			 35,   0,   0,   0,   4,   0,
			  0,   1,   4,   3,   0,   0,
			  0,   2,   2,   1,   1,   1
    	};
		
		//the tensor describes layers, rows, columns, source of input
		utility::fixedsizearray1< Size, 3 > dims;
		dims[1] = 5; dims[2] = 6; dims[3] = 6;
		MathNTensor< Real, 3 > const input_tensor( dims, values);
		MathTensor< Real > const input_tensor2( 5, 6, 6, values);
		
		BorderFlag border2[3] = { e_Natural, e_Natural, e_Periodic};
		
		//these vectors are used to input the starting points and the
		//grid width delta of every dimension (x, y, z) into the spline
		const double start2[] = {10, 10, 10};
		const double delta2[] = {10, 10, 10};
		
		const bool lin_cont2[3] = { true, true, true};
		
		const std::pair< double, double> first_be2[3] =
		{
			std::pair< double, double>( 10, 10),
			std::pair< double, double>( 10, 10),
			std::pair< double, double>( 10, 10)
		};
		
		utility::fixedsizearray1< BorderFlag, 3 > border( e_Natural );
		border[ 3 ] = e_Periodic;
		utility::fixedsizearray1< double, 3> start( 10 );
		utility::fixedsizearray1< double, 3> delta( 10 );
		utility::fixedsizearray1< bool, 3> lin_cont( true );
		utility::fixedsizearray1< std::pair< double, double>, 3> first_be( std::pair< double, double>( 10, 10) );
		
		
		
		
		// natural means, the second order derivative at beginning and end is 0
		// this leads to a rather slow changing first order derivative meaning
		// "nearly linear" behavior of the spline
		// in z-direction the spline is made periodic
		PolycubicSpline<3> naturalspline;
		TricubicSpline naturalspline2;
		naturalspline.train( border, start, delta, input_tensor, lin_cont, first_be);
		naturalspline2.train( border2, start2, delta2, input_tensor2, lin_cont2, first_be2);
		
		for( double x( 19.9); x < 20.1; x += 0.01)
		{
			utility::fixedsizearray1< Real, 3 > vals( 10.0 );
			vals[1] = x;
			//std::cout << x << " " << naturalspline.F( vals ) << " " << naturalspline.dFdxi( 1, vals ) << std::endl;
			//std::cout << x << " " << naturalspline2.F( x, 10, 10) << " " << naturalspline2.dFdx( x, 10, 10) << std::endl;
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline2.F( x, 10, 10 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 1, vals ), naturalspline2.dFdx( x, 10, 10 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 2, vals ), naturalspline2.dFdy( x, 10, 10 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 3, vals ), naturalspline2.dFdz( x, 10, 10 ), 1e-12 );
		}
		
		for( double x( 49); x < 50.5; x += 0.1)
		{
			utility::fixedsizearray1< Real, 3 > vals( 10.0 );
			vals[1] = x;
			
			std::cout << x << " ";
			std::cout << naturalspline.F( vals ) << " ";
			std::cout << naturalspline.dFdxi( 1, vals ) << std::endl;
			
			std::cout << x << " ";
			std::cout << naturalspline2.F( x, 10, 10 ) << " ";
			std::cout << naturalspline2.dFdx( x, 10, 10 ) << std::endl;
			
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline2.F( x, 10, 10 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 1, vals ), naturalspline2.dFdx( x, 10, 10 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 2, vals ), naturalspline2.dFdy( x, 10, 10 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 3, vals ), naturalspline2.dFdz( x, 10, 10 ), 1e-12 );
		}
		
		for( double x( 49); x < 50.5; x += 0.1)
		{
			utility::fixedsizearray1< Real, 3 > vals( 15.0 );
			vals[1] = x;
			
			std::cout << x << " ";
			std::cout << naturalspline.F( vals ) << " ";
			std::cout << naturalspline.dFdxi( 1, vals ) << std::endl;
			
			std::cout << x << " ";
			std::cout << naturalspline2.F( x, 15, 15 ) << " ";
			std::cout << naturalspline2.dFdx( x, 15, 15 ) << std::endl;
			
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline2.F( x, 15, 15 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 1, vals ), naturalspline2.dFdx( x, 15, 15 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 2, vals ), naturalspline2.dFdy( x, 15, 15 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 3, vals ), naturalspline2.dFdz( x, 15, 15 ), 1e-12 );
		}
		
		std::cout << "To show periodic behavior in z-direction" << std::endl;
		std::cout << " z        F(10, 10, z)       F(10, 10, z+60)" << std::endl;
		for( double z( 20); z < 60; z += 4)
		{
			utility::fixedsizearray1< Real, 3 > vals( 10.0 );
			utility::fixedsizearray1< Real, 3 > valsp60( 10.0 );
			vals[3] = z;
			valsp60[3] = z+60;
			
			std::cout << z << " ";
			std::cout << naturalspline.F( vals ) << " ";
			std::cout <<  naturalspline.F( valsp60 ) << std::endl;
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline.F( valsp60 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 1, vals ), naturalspline.dFdxi( 1, valsp60 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 2, vals ), naturalspline.dFdxi( 2, valsp60 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 3, vals ), naturalspline.dFdxi( 3, valsp60 ), 1e-12 );
			std::cout << z << " ";
			std::cout << naturalspline2.F( 10, 10, z) << " ";
			std::cout <<  naturalspline2.F( 10, 10, z+60) << std::endl;
			TS_ASSERT_DELTA( naturalspline2.F( 10, 10, z), naturalspline2.F( 10, 10, z+60), 1e-12 );
			TS_ASSERT_DELTA( naturalspline2.dFdx( 10, 10, z ), naturalspline2.dFdx( 10, 10, z+60 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline2.dFdy( 10, 10, z ), naturalspline2.dFdy( 10, 10, z+60 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline2.dFdz( 10, 10, z ), naturalspline2.dFdz( 10, 10, z+60 ), 1e-12 );
			
			
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline2.F( 10, 10, z+60), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 1, vals ), naturalspline2.dFdx( 10, 10, z+60 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 2, vals ), naturalspline2.dFdy( 10, 10, z+60 ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 3, vals ), naturalspline2.dFdz( 10, 10, z+60 ), 1e-12 );
		}
		
		std::cout << "space diagonal values: " << std::endl;
		for( double r( -10); r < 100; r += 0.1)
		{
			utility::fixedsizearray1< Real, 3 > vals( r );
			
			std::cout << r << " ";
			std::cout << naturalspline.F( vals ) << " ";
			std::cout << naturalspline2.F( r, r, r) << " ";
			
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline2.F( r, r, r ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 1, vals ), naturalspline2.dFdx( r, r, r ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 2, vals ), naturalspline2.dFdy( r, r, r ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 3, vals ), naturalspline2.dFdz( r, r, r ), 1e-12 );
		}
		
		
		std::cout << "corroborate around 46 20 10" << std::endl;
		for( double z( 9); z < 11; z += 0.1)
		{
			utility::fixedsizearray1< Real, 3 > vals( 10.0 );
			vals[1] = 46;
			vals[2] = 20;
			vals[3] = z;
			
			TS_ASSERT_DELTA( naturalspline.F( vals ), naturalspline2.F( 46, 20, z), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 1, vals ), naturalspline2.dFdx( 46, 20, z ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 2, vals ), naturalspline2.dFdy( 46, 20, z ), 1e-12 );
			TS_ASSERT_DELTA( naturalspline.dFdxi( 3, vals ), naturalspline2.dFdz( 46, 20, z ), 1e-12 );
		}
		
		//this example describes a function f(x, y, z), where (x, y, z) resembles a helix
		//the spline is periodic in every direction, but this is not shown here
		PolycubicSpline< 3 > ps;
		TricubicSpline ts;
		
		//std::cout <<  "Training new purely periodic spline" << std::endl;
		
		border[ 1 ] = e_Periodic;
		border[ 2 ] = e_Periodic;
		border[ 3 ] = e_Periodic;
		
		border2[ 0 ] = e_Periodic;
		border2[ 1 ] = e_Periodic;
		border2[ 2 ] = e_Periodic;


		ts.train( border2, start2, delta2, input_tensor2, lin_cont2, first_be2);
		ps.train( border, start, delta, input_tensor, lin_cont, first_be);
		
		//last function example
		//std::cout <<  "Function values along a helix" << std::endl;
		//std::cout << "x F( 10+x, 10+10*std::cos( pi / 18 * x), 10+10*std::sin( pi / 18 * x))" << std::endl;
		for( double x( 0); x <= 36; ++x)
		{
			utility::fixedsizearray1< Real, 3 > vals;
			vals[1] = 10 + x;
			vals[2] = 10 + 10 * std::cos( numeric::constants::d::pi / 18 * x);
			vals[3] = 10 + 10 * std::sin( numeric::constants::d::pi / 18 * x);
			
			TS_ASSERT_DELTA(  ps.F( vals ), ts.F( 10 + x, 10 + 10 * std::cos( numeric::constants::d::pi / 18 * x), 10 + 10 * std::sin( numeric::constants::d::pi / 18 * x)), 1e-12 );
			TS_ASSERT_DELTA(  ps.dFdxi( 1, vals ), ts.dFdx( 10 + x, 10 + 10 * std::cos( numeric::constants::d::pi / 18 * x), 10 + 10 * std::sin( numeric::constants::d::pi / 18 * x)), 1e-12 );
			TS_ASSERT_DELTA(  ps.dFdxi( 2, vals ), ts.dFdy( 10 + x, 10 + 10 * std::cos( numeric::constants::d::pi / 18 * x), 10 + 10 * std::sin( numeric::constants::d::pi / 18 * x)), 1e-12 );
			TS_ASSERT_DELTA(  ps.dFdxi( 3, vals ), ts.dFdz( 10 + x, 10 + 10 * std::cos( numeric::constants::d::pi / 18 * x), 10 + 10 * std::sin( numeric::constants::d::pi / 18 * x)), 1e-12 );

		}
	}

};

