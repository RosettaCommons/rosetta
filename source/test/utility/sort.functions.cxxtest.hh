// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/sort.functions.cxxtest.hh
/// @brief  sort.functions.cxxtest: test suite for utility::sort.functions
/// @author Brian Coventry (bcov@uw.edu)

// Testing headers
#include <cxxtest/TestSuite.h>


// Project headers
#include <utility/sort.functions.hh>
#include <core/types.hh>

// C/C++ headers
#include <iostream>
#include <string>
#include <sstream>

//#include <basic/Tracer.hh>

// static basic::Tracer TR( "utility.sort.functions.cxxtest.hh" );

class SortFunctionsTest : public CxxTest::TestSuite {
public:

	void test_argsort() {

		utility::vector1< core::Size > to_sort { 0, 2, 5, 7, 1, 4 };
		utility::vector1< core::Size > sorted = utility::argsort( to_sort );

		TS_ASSERT_EQUALS( sorted[1], 1 );
		TS_ASSERT_EQUALS( sorted[2], 5 );
		TS_ASSERT_EQUALS( sorted[3], 2 );
		TS_ASSERT_EQUALS( sorted[4], 6 );
		TS_ASSERT_EQUALS( sorted[5], 3 );
		TS_ASSERT_EQUALS( sorted[6], 4 );
	}

	void test_fractional_rank() {

		{
			utility::vector1< core::Size > to_sort { 0, 2, 5, 7, 1, 4 };
			utility::vector1< core::Real > sorted = utility::fractional_rank( to_sort );

			TS_ASSERT_DELTA( sorted[1], 1, 0.001 );
			TS_ASSERT_DELTA( sorted[2], 3, 0.001 );
			TS_ASSERT_DELTA( sorted[3], 5, 0.001 );
			TS_ASSERT_DELTA( sorted[4], 6, 0.001 );
			TS_ASSERT_DELTA( sorted[5], 2, 0.001 );
			TS_ASSERT_DELTA( sorted[6], 4, 0.001 );
		}

		{
			utility::vector1< core::Size > to_sort { 0, 4, 2, 4, 2, 4 };
			utility::vector1< core::Real > sorted = utility::fractional_rank( to_sort );

			TS_ASSERT_DELTA( sorted[1], 1, 0.001 );
			TS_ASSERT_DELTA( sorted[2], 5, 0.001 );
			TS_ASSERT_DELTA( sorted[3], 2.5, 0.001 );
			TS_ASSERT_DELTA( sorted[4], 5, 0.001 );
			TS_ASSERT_DELTA( sorted[5], 2.5, 0.001 );
			TS_ASSERT_DELTA( sorted[6], 5, 0.001 );
		}

	}

};
