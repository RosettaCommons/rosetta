// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/Quaterion.cxxtest.hh
/// @brief  test suite for numeric::Quaternion
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Core headers
#include <core/types.hh>

// Unit headers
#include <numeric/interpolation/Histogram.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <fstream>


// --------------- Test Class --------------- //

using namespace core;
using utility::vector1;

class HistogramTests : public CxxTest::TestSuite {

typedef numeric::interpolation::Histogram<Real, Real> Histogram;

public:

	// Shared data elements
	vector1<Real> data;

	HistogramTests():
		data() {}

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		Real values[] = {1.0, 0.0, 2.0, 2.0, -3.5}; //Match contents of Histogram_sample.hist
		data.assign(values, values+5);
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	/// @brief test construction
	void test_Histogram_Constructors() {
		Histogram h1(data, -5.0, 2.0);

		//Don't test the file input. Including gives weird linker errors.
		std::ifstream sample_hist("numeric/interpolation/Histogram_sample.hist");
		TS_ASSERT(sample_hist.good());
		Histogram h2(sample_hist);

		Histogram h3(h1);

		TS_ASSERT_EQUALS(data,h1.densities());
		TS_ASSERT_EQUALS(data,h2.densities());
		TS_ASSERT_EQUALS(data,h3.densities());

		TS_ASSERT_EQUALS(h2.minimum(),-5.0);
		TS_ASSERT_EQUALS(h2.maximum(),5.0);
		TS_ASSERT_EQUALS(h2.step_size(),2.0);
		TS_ASSERT(! h2.periodic() );
		TS_ASSERT_EQUALS(h2.bin_placement(), Histogram::left);
		TS_ASSERT_EQUALS(h2.interpolator(), Histogram::flat);

	}

	/// @brief test that the bounds are correct for the various BinPlacements
	void test_Histogram_bounds_flat() {
		Histogram h_left(data, -5.0, 2.0, false, Histogram::left, Histogram::flat);

		TS_ASSERT_EQUALS(h_left.minimum(), -5.0);
		TS_ASSERT_EQUALS(h_left.maximum(), 5.0);
		TS_ASSERT_EQUALS(h_left.first_bin(), -5.0);
		TS_ASSERT_EQUALS(h_left.last_bin(), 3.0);

		Histogram h_center(data, -5.0, 2.0, false, Histogram::center,Histogram::flat);

		TS_ASSERT_EQUALS(h_center.minimum(), -5.0);
		TS_ASSERT_EQUALS(h_center.maximum(), 5.0);
		TS_ASSERT_EQUALS(h_center.first_bin(), -5.0);
		TS_ASSERT_EQUALS(h_center.last_bin(), 3.0);
	}
	/// @brief test that the bounds are correct for the various BinPlacements
	void test_Histogram_bounds_linear() {
		Histogram h_left(data, -5.0, 2.0, false, Histogram::left, Histogram::linear);

		TS_ASSERT_EQUALS(h_left.minimum(), -5.0);
		TS_ASSERT_EQUALS(h_left.maximum(), 3.0);
		TS_ASSERT_EQUALS(h_left.first_bin(), -5.0);
		TS_ASSERT_EQUALS(h_left.last_bin(), 3.0);

		Histogram h_center(data, -5.0, 2.0, false, Histogram::center,Histogram::linear);

		TS_ASSERT_EQUALS(h_center.minimum(), -4.0);
		TS_ASSERT_EQUALS(h_center.maximum(), 4.0);
		TS_ASSERT_EQUALS(h_center.first_bin(), -5.0);
		TS_ASSERT_EQUALS(h_center.last_bin(), 3.0);
	}

	void test_Histogram_interpolation_flat() {

		bool in_range;
		Real y = 0.0 ;

		// Left
		Histogram h_left(data, -5.0, 2.0, false, Histogram::left, Histogram::flat);

		// Valid range
		in_range = h_left.interpolate(-4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_left.interpolate(-1.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 2.0);

		//Lower boundary
		in_range = h_left.interpolate(-5.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_left.interpolate(-5.1, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_left.interpolate(-100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);

		//Upper boundary
		in_range = h_left.interpolate(4.9, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -3.5);
		in_range = h_left.interpolate(5.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);
		in_range = h_left.interpolate(100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);


		// Center
		Histogram h_center(data, -5.0, 2.0, false, Histogram::center, Histogram::flat);

		// Valid range
		in_range = h_center.interpolate(-4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_center.interpolate(-1.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 2.0);

		//Lower boundary
		in_range = h_center.interpolate(-5.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_center.interpolate(-5.1, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_center.interpolate(-100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);

		//Upper boundary
		in_range = h_center.interpolate(4.9, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -3.5);
		in_range = h_center.interpolate(5.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);
		in_range = h_center.interpolate(100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);

	}

	void test_Histogram_interpolation_linear() {

		bool in_range;
		Real y = 0.0 ;

		// Left
		Histogram h_left(data, -5.0, 2.0, false, Histogram::left, Histogram::linear);

		// Valid range
		in_range = h_left.interpolate(-4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.5);
		in_range = h_left.interpolate(0.75, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 2.0);
		in_range = h_left.interpolate(-4.5, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.75);
		in_range = h_left.interpolate(1+4.0/5.5, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, 0.0, 1e-8);

		// Lower boundary
		in_range = h_left.interpolate(-5.01, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_left.interpolate(-100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);

		// Upper boundary
		in_range = h_left.interpolate(3.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);
		in_range = h_left.interpolate(100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);

		// Center
		Histogram h_center(data, -6.0, 2.0, false, Histogram::center, Histogram::linear);

		// Valid range
		in_range = h_center.interpolate(-5.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_center.interpolate(-4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.5);
		in_range = h_center.interpolate(0.75, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 2.0);
		in_range = h_center.interpolate(-4.5, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.75);
		in_range = h_center.interpolate(1+4.0/5.5, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, 0.0,1e-8);

		// Lower boundary
		in_range = h_center.interpolate(-5.01, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_center.interpolate(-100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 1.0);

		// Upper boundary
		in_range = h_center.interpolate(3.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);
		in_range = h_center.interpolate(100.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, -3.5);
	}

	void test_Histogram_interpolation_periodic_linear() {


		bool in_range;
		Real y = 0.0 ;

		// Left
		Histogram h_left(data, -5.0, 2.0, true, Histogram::left, Histogram::linear);

		TS_ASSERT_EQUALS(h_left.minimum(),-5.0);
		TS_ASSERT_EQUALS(h_left.maximum(),3.0);
		TS_ASSERT_EQUALS(h_left.first_bin(),-5.0);
		TS_ASSERT_EQUALS(h_left.last_bin(),3.0);

		// Valid range. Period 10, so adding multiples of 10 does nothing
		in_range = h_left.interpolate(-4.0+10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.5);
		in_range = h_left.interpolate(0.75+100, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 2.0);
		in_range = h_left.interpolate(-4.5-50, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.75);
		in_range = h_left.interpolate(1+4.0/5.5-10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, 0.0,1e-8);

		// Boundary
		in_range = h_left.interpolate(4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, -3.5+4.5/2.0,1e-8);
		in_range = h_left.interpolate(5.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_left.interpolate(-6.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, -3.5+4.5/2.0,1e-8);

		// Center
		Histogram h_center(data, -6.0, 2.0, true, Histogram::center, Histogram::linear);

		TS_ASSERT_EQUALS(h_center.minimum(),-5.0);
		TS_ASSERT_EQUALS(h_center.maximum(),3.0);
		TS_ASSERT_EQUALS(h_center.first_bin(),-6.0);
		TS_ASSERT_EQUALS(h_center.last_bin(),2.0);

		// Valid range. Period 10, so adding multiples of 10 does nothing
		in_range = h_center.interpolate(-4.0+10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.5);
		in_range = h_center.interpolate(0.75+100, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 2.0);
		in_range = h_center.interpolate(-4.5-50, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.75);
		in_range = h_center.interpolate(1+4.0/5.5-10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, 0.0,1e-8);
		in_range = h_center.interpolate(3.0-30, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -3.5);

		// Boundary
		in_range = h_center.interpolate(4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, -3.5+4.5/2.0,1e-8);
		in_range = h_center.interpolate(5.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 1.0);
		in_range = h_center.interpolate(-6.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_DELTA(y, -3.5+4.5/2.0,1e-8);
	}

	void test_derivative_linear() {
		bool in_range;
		Real y = 0.0 ;

		// Left
		Histogram h_left(data, -5.0, 2.0, false, Histogram::left, Histogram::linear);

		in_range = h_left.derivative(-5.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_left.derivative(-4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_left.derivative(0.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.0);
		in_range = h_left.derivative(2.9, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -5.5/2.0);

		// Bounds
		in_range = h_left.derivative(-5.1, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 0.0);
		in_range = h_left.derivative(3.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 0.0);

		// Center
		Histogram h_center(data, -6.0, 2.0, false, Histogram::center, Histogram::linear);

		in_range = h_center.derivative(-5.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_center.derivative(-4.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_center.derivative(0.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.0);
		in_range = h_center.derivative(2.9, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -5.5/2.0);

		// Bounds
		in_range = h_center.derivative(-5.1, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 0.0);
		in_range = h_center.derivative(3.0, y);
		TS_ASSERT(!in_range);
		TS_ASSERT_EQUALS(y, 0.0);
	}

	void test_derivative_periodic_linear() {
		bool in_range;
		Real y = 0.0 ;

		// Left
		Histogram h_left(data, -5.0, 2.0, true, Histogram::left, Histogram::linear);

		in_range = h_left.derivative(-5.0+100, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_left.derivative(-4.0-20, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_left.derivative(0.0-10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.0);
		in_range = h_left.derivative(2.9+10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -5.5/2.0);

		// Bounds
		in_range = h_left.derivative(3.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 4.5/2.0);
		in_range = h_left.derivative(-6.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 4.5/2.0);
		in_range = h_left.derivative(3.41+40, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 4.5/2.0);

		// center
		Histogram h_center(data, -6.0, 2.0, true, Histogram::center, Histogram::linear);

		in_range = h_center.derivative(-5.0+100, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_center.derivative(-4.0-20, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -.5);
		in_range = h_center.derivative(0.0-10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 0.0);
		in_range = h_center.derivative(2.9+10, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, -5.5/2.0);

		// Bounds
		in_range = h_center.derivative(3.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 4.5/2.0);
		in_range = h_center.derivative(-6.0, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 4.5/2.0);
		in_range = h_center.derivative(3.41+40, y);
		TS_ASSERT(in_range);
		TS_ASSERT_EQUALS(y, 4.5/2.0);
	}


};

