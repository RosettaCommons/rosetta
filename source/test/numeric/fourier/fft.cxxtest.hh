// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/Quaterion.cxxtest.hh
/// @brief  test suite for numeric::Quaternion
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/fourier/FFT.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <complex>

class FFTTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	/// @brief Constructors
	void test_small_transform() {
		// small 1d r->c transform versus stored result
		ObjexxFCL::FArray1D< double > x;
		ObjexxFCL::FArray1D< std::complex< double > > fx;
		x.dimension(6);
		for (int i=0; i<6; ++i) {
			x[i] = i;
		}
		numeric::fourier::fft ( x, fx );
		TS_ASSERT_DELTA( fx[0].real(), 15.0000000, 1e-6);
		TS_ASSERT_DELTA( fx[0].imag(),  0.0000000, 1e-6);
		TS_ASSERT_DELTA( fx[1].real(), -3.0000000, 1e-6);
		TS_ASSERT_DELTA( fx[1].imag(),  5.1961524, 1e-6);
		TS_ASSERT_DELTA( fx[2].real(), -3.0000000, 1e-6);
		TS_ASSERT_DELTA( fx[2].imag(),  1.7320508, 1e-6);
		TS_ASSERT_DELTA( fx[3].real(), -3.0000000, 1e-6);
		TS_ASSERT_DELTA( fx[3].imag(),  0.0000000, 1e-6);
	}

	void test_1D_transform() {
		ObjexxFCL::FArray1D< std::complex< double > > x, y;
		ObjexxFCL::FArray1D< std::complex< double > > fx;

		// dim==120 tests bfly2,3,4,5
		x.dimension(120);
		for (int i=0; i<120; ++i) {
			x[i] = std::complex< double >( 2*i, 2*i+1);
		}
		numeric::fourier::fft ( x, fx );
		numeric::fourier::ifft( fx, y );

		double err = 0.0;
		for (int i=0; i<120; ++i)
			err += std::norm(y[i]-x[i]);
		TS_ASSERT_DELTA( err, 0.0, 1e-6);
	}

	/// @brief test normalization
	void test_2D_transform() {
		ObjexxFCL::FArray2D< std::complex< double > > x, y;
		ObjexxFCL::FArray2D< std::complex< double > > fx;

		// dim==120 tests bfly2,3,4,5
		x.dimension(120,120);
		for (int i=0; i<120*120; ++i) {
			x[i] = std::complex< double >( 2*i, 2*i+1);
		}
		numeric::fourier::fft2 ( x, fx );
		numeric::fourier::ifft2( fx, y );

		double err = 0.0;
		for (int i=0; i<120*120; ++i)
			err += std::norm(y[i]-x[i]);
		TS_ASSERT_DELTA( err, 0.0, 1e-6);
	}

	/// @brief test conjugates
	void test_3D_transform() {
		ObjexxFCL::FArray3D< std::complex< double > > x, y;
		ObjexxFCL::FArray3D< std::complex< double > > fx;

		// dim==120 tests bfly2,3,4,5
		x.dimension(120,120,120);
		for (int i=0; i<120*120*120; ++i) {
			x[i] = std::complex< double >( 2*i, 2*i+1);
		}
		numeric::fourier::fft3 ( x, fx );
		numeric::fourier::ifft3( fx, y );

		double err = 0.0;
		for (int i=0; i<120*120*120; ++i)
			err += std::norm(y[i]-x[i]);
		TS_ASSERT_DELTA( err, 0.0, 1e-6);
	}

};


