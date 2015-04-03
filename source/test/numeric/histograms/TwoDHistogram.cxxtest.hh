// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Test headers
#include <cxxtest/TestSuite.h>


#include <numeric/histograms/TwoDHistogram.hh>

// --------------- Test Class --------------- //
class TwoDHistogramTests : public CxxTest::TestSuite{

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

	} //Match contents of Histogram_sample.hist


	// Shared finalization goes here.
	void tearDown() {

	}


	void test_two_d_histogram(){
		numeric::histograms::TwoDHistogram<platform::Real, platform::Real> histogram_real;

		numeric::histograms::TwoDHistogram<platform::Size, platform::Size> histogram_size;


		//make key
		std::pair<platform::Size, platform::Size> key_values(1,2);
		std::pair<platform::Real, platform::Real> key_values_real(1.5,2.5);

		//make counts
		platform::Size counts(50);


		histogram_size.insert_data(1,2, counts);
		histogram_real.insert_data(1.5,2.5, counts);


		TS_ASSERT_EQUALS(counts, histogram_size.lookup_counts(1,2));
		TS_ASSERT_EQUALS(counts, histogram_real.lookup_counts(1.5,2.5));


	}


};
