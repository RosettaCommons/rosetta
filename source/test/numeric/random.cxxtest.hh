// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random.cxxtest.hh
/// @brief  test suite for numeric::random
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/random/random.hh>

#include <vector>
#include <iostream>
#include <sstream>

#include <test/UTracer.hh>


class RandomSystemTests : public CxxTest::TestSuite
{
public:
	RandomSystemTests() {}

	// Shared initialization goes here.
	void setUp() {
		numeric::random::rg().set_seed( "mt19937", 999 );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// ------------------------------------------ //
	/// @brief test how Uniform distribution is
	void test_Uniform() {
		numeric::random::rg().set_seed("mt19937", 999);

		const int size = 5000;  ///< size of array where we store all uniform numbers
		const int N = 20100100; //< number of trials
		std::vector<int> V(size);
		for(unsigned int i=0; i<V.size(); i++) V[i] = 0;

		int min_ind = size/2;
		int max_ind = size/2;
		for(int i=0; i<N; i++) {
			double r = numeric::random::rg().uniform();
			int ind = int( r * double(size-1e-10) );
			if( min_ind > ind ) {
				min_ind = ind;
				//std::cout << ind << " ";
			}
			if( max_ind < ind ) {
				max_ind = ind;
				//std::cout << ind << " ";
			}
			V[ind] += 1;
		}
		TS_ASSERT_EQUALS( min_ind, 0 );
		TS_ASSERT_EQUALS( max_ind, size-1 );

		int min_c = N;
		int max_c = 0;
		for(unsigned int i=0; i<V.size(); i++) {
			if( min_c > V[i] ) min_c = V[i];
			if( max_c < V[i] ) max_c = V[i];
		}
		double rate = double(max_c)/min_c - 1.;
		//std::cout << min_c << " ";
		//std::cout << max_c << " ";
		//std::cout << rate << " ";
		TS_ASSERT_LESS_THAN(rate, .15);
		//TS_TRACE("some message");
	}

	void test_SaveAndRestore() {
		const int size = 5000;
		std::vector<double> V1(size), V2(size);

		std::ostringstream oss;
		numeric::random::rg().saveState(oss);
		for(int i = 0; i < size; i++) {
			V1[i] = numeric::random::rg().uniform();
		}

		std::istringstream iss( oss.str() );
		numeric::random::rg().restoreState(iss);
		for(int i = 0; i < size; i++) {
			V2[i] = numeric::random::rg().uniform();
		}

		TS_ASSERT_EQUALS(V1, V2); // vector == is element-wise
	}

	// -------------------------------------------
	/// @brief test that number from generator mt19937 are the same
	void test_mt19937() {
		const int NumbersToTest = 10000;
		numeric::random::rg().set_seed("mt19937", 999);

		// UTracer log file
		test::UTracer  UT("numeric/mt19937.u");

		UT << std::setprecision(16);

		for(int i=0; i<NumbersToTest; i++) {
			UT << numeric::random::rg().uniform() << std::endl;
		}
	}
};

