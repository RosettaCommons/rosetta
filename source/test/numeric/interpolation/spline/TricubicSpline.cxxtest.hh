// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <numeric/interpolation/spline/TricubicSpline.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>
#include <numeric/MathVector_operations.hh>



// --------------- Test Class --------------- //


class TricubicSplineTests : public CxxTest::TestSuite {



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

	void test_tricubic_spline_data_access() {

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
		
		//the tensor describes layers, rows, columns, source of input
		const MathTensor< Real > input_tensor( 5, 6, 6, values);
		
		BorderFlag border[3] = { e_Natural, e_Natural, e_Periodic};
		
		//these vectors are used to input the starting points and the
		//grid width delta of every dimension (x, y, z) into the spline
		const double start[] = {10, 10, 10};
		const double delta[] = {10, 10, 10};
		
		const bool lin_cont[3] = { true, true, true};
		
		//this vector controls the behavior of the spline at the beginning and
		//end of every dimension, only has impact for BorderFlag FIRSTDER
		
		//every pair describes the value of the first order derivative at
		//start and end
		const std::pair< double, double> first_be[3] =
		{
			std::pair< double, double>( 10, 10),
			std::pair< double, double>( 10, 10),
			std::pair< double, double>( 10, 10)
		};
		
		
		// natural means, the second order derivative at beginning and end is 0
		// this leads to a rather slow changing first order derivative meaning
		// "nearly linear" behavior of the spline
		// in z-direction the spline is made periodic
		TricubicSpline naturalspline;
		naturalspline.train( border, start, delta, input_tensor, lin_cont, first_be);
		
		//the next line generates a class from a trained spline that can give back F(x), dF(x) and FdF(x)
		//naturalspline.WriteCodeInFile("natural_trained", GetExampleSourceCodePath());
		
		//BCL_Message
		//(
		 //util::Message::e_Standard,
		 //std::cout <<
		 //"Example for a spline with natural behavior in x- and y-direction (F_xx(x_0, y, z)=F_xx(x_dim-1, y, z)= "
		 //"F_yy(x, y_0, z)=F_yy(x, y_dim-1, z)=0) "
		 //"and periodic behavior in z-direction (F(x, y, z + n * dimz * deltaz)=F(x, y, z)) :" << std::endl;
		
		//BCL_Message( util::Message::e_Standard,
		//std::cout << "To show continuous behavior at the end of a cell" << std::endl;
		//BCL_Message( util::Message::e_Standard, " x        F(x, 10, 10)       F_x(x, 10, 10)");
		//std::cout << " x        F(x, 10, 10)       F_x(x, 10, 10)" << std::endl;
		for( double x( 19.9); x < 20.1; x += 0.01)
		{
			//std::cout << x << " " << naturalspline.F( x, 10, 10) << " " << naturalspline.dFdx( x, 10, 10) << std::endl;
			//BCL_Message( util::Message::e_Standard, fmt2( x) + fmt( naturalspline.F( x, 10, 10)) + fmt( naturalspline.dFdx( x, 10, 10)));
			//BCL_Message
			//(
			//util::Message::e_Standard, fmt2(x) + fmt( model::TricubicSpline_natural_trained().F
			//														( math::MakeVector( x, 10., 10.)))
			// + fmt( model::TricubicSpline_natural_trained().dFdx( math::MakeVector( x, 10., 10.)))
			// );
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
			//std::cout <<  z << " ";
			//std::cout << naturalspline.F( 10, 10, z) << " ";
			//std::cout <<  naturalspline.F( 10, 10, z+60) << std::endl;
			TS_ASSERT_DELTA( naturalspline.F( 10, 10, z), naturalspline.F( 10, 10, z+60), 1e-12 );
		}
		
		//this example describes a function f(x, y, z), where (x, y, z) resembles a helix
		//the spline is periodic in every direction, but this is not shown here
		TricubicSpline ts;
		
		//std::cout <<  "Training new purely periodic spline" << std::endl;
		
		border[ 0 ] = e_Periodic;
		border[ 1 ] = e_Periodic;
		
		ts.train( border, start, delta, input_tensor, lin_cont, first_be);
		
		//last function example
		//std::cout <<  "Function values along a helix" << std::endl;
		//std::cout << "x F( 10+x, 10+10*std::cos( pi / 18 * x), 10+10*std::sin( pi / 18 * x))" << std::endl;
		for( double x( 0); x <= 36; ++x)
		{
			//std::cout <<  x << " ";
			//std::cout << ts.F( 10 + x, 10 + 10 * std::cos( numeric::constants::d::pi / 18 * x), 10 + 10 * std::sin( numeric::constants::d::pi / 18 * x)) << std::endl;
		}
		
	}


};
