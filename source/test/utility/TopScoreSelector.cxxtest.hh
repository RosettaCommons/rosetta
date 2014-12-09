// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/TopScoreSelector.cxxtest.hh
/// @brief  test suite for utility::TopScoreSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Package headers
#include <utility/TopScoreSelector.hh>
#include <cxxtest/TestSuite.h>

// Project headers
#include <core/types.hh>

/// C++ headers
#include <iostream>


class TopScoreSelectorTests : public CxxTest::TestSuite {

private:

public:

	void setUp() {
	}

	/// @brief Constructor Test
	void test_TopScoreSelector_constructor() {

		utility::TopScoreSelector< int, core::Real > ts_selector;
		ts_selector.n_to_keep( 5 );
		TS_ASSERT( ts_selector.n_to_keep() == 5 );
		TS_ASSERT( ts_selector.low_is_better() );
		TS_ASSERT( ts_selector.better( 2.5, 5.5 ) );
	}

	void test_TopScoreSelector_insert_one() {

		utility::TopScoreSelector< int, core::Real > ts_selector;
		ts_selector.n_to_keep( 5 );
		ts_selector.insert( 4, 19.0 );

		TS_ASSERT( ts_selector.size() == 1 );
		TS_ASSERT( ts_selector[ 1 ] == 4 );

	}

	void test_TopScoreSelector_insert_five() {

		utility::TopScoreSelector< int, core::Real > ts_selector;
		ts_selector.n_to_keep( 5 );

		ts_selector.insert( 4, 19.0 );
		ts_selector.insert( 3, 10.5 );
		ts_selector.insert( 12, 20.0 );
		ts_selector.insert( 50, 13.0 );
		ts_selector.insert( -1, 15.0 );

		TS_ASSERT( ts_selector.size() == 5 );
		TS_ASSERT( ts_selector[ 1 ] == 3 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == 50 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == -1 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 4 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == 12 ); //std::cout << ts_selector[ 5 ] << std::endl;

	}

	void test_TopScoreSelector_insert_fiveB() {

		utility::TopScoreSelector< int, core::Real > ts_selector;
		ts_selector.n_to_keep( 5 );

		ts_selector.insert( 4, 19.0 );
		ts_selector.insert( 3, 10.5 );
		ts_selector.insert( 12, 20.0 );
		ts_selector.insert( 50, 13.0 );
		ts_selector.insert( -1, 5.0 );

		TS_ASSERT( ts_selector.size() == 5 );
		TS_ASSERT( ts_selector[ 1 ] == -1 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == 3 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == 50 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 4 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == 12 ); //std::cout << ts_selector[ 5 ] << std::endl;

	}

	void test_TopScoreSelector_insert_sixB() {

		utility::TopScoreSelector< int, core::Real > ts_selector;
		ts_selector.n_to_keep( 5 );

		ts_selector.insert( 4, 19.0 );
		ts_selector.insert( 3, 10.5 );
		ts_selector.insert( 12, 20.0 );
		ts_selector.insert( 50, 13.0 );
		ts_selector.insert( -1, 5.0 );

		TS_ASSERT( ts_selector.size() == 5 );
		TS_ASSERT( ts_selector[ 1 ] == -1 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == 3 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == 50 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 4 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == 12 ); //std::cout << ts_selector[ 5 ] << std::endl;

		ts_selector.insert( 66, -3.0 );
		TS_ASSERT( ts_selector[ 1 ] == 66 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == -1 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == 3 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 50 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == 4 ); //std::cout << ts_selector[ 5 ] << std::endl;

	}

	void test_TopScoreSelector_replace_initial_five() {

		utility::TopScoreSelector< int, core::Real > ts_selector;
		ts_selector.n_to_keep( 5 );

		ts_selector.insert( 4, 19.0 );
		ts_selector.insert( 3, 10.5 );
		ts_selector.insert( 12, 20.0 );
		ts_selector.insert( 50, 13.0 );
		ts_selector.insert( -1, 5.0 );

		TS_ASSERT( ts_selector.size() == 5 );
		TS_ASSERT( ts_selector[ 1 ] == -1 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == 3 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == 50 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 4 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == 12 ); //std::cout << ts_selector[ 5 ] << std::endl;

		ts_selector.insert( 104, -19.0 );
		ts_selector.insert( 103, -10.5 );
		ts_selector.insert( 112, -20.0 );
		ts_selector.insert( 150, -13.0 );
		ts_selector.insert( -101, -5.0 );

		TS_ASSERT( ts_selector[ 1 ] == 112 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == 104 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == 150 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 103 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == -101 ); //std::cout << ts_selector[ 5 ] << std::endl;

	}

	void test_TopScoreSelector_replace_initial_five_w_higher_being_better() {

		utility::TopScoreSelector< int, core::Real > ts_selector;
		ts_selector.n_to_keep( 5 );
		ts_selector.low_is_better( false );

		ts_selector.insert( 4, 19.0 );
		ts_selector.insert( 3, 10.5 );
		ts_selector.insert( 12, 20.0 );
		ts_selector.insert( 50, 13.0 );
		ts_selector.insert( -1, 5.0 );

		TS_ASSERT( ts_selector.size() == 5 );
		TS_ASSERT( ts_selector[ 1 ] == -1 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == 3 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == 50 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 4 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == 12 ); //std::cout << ts_selector[ 5 ] << std::endl;

		ts_selector.insert( 104, -19.0 );
		ts_selector.insert( 103, -10.5 );
		ts_selector.insert( 112, -20.0 );
		ts_selector.insert( 150, -13.0 );
		ts_selector.insert( -101, -5.0 );

		TS_ASSERT( ts_selector[ 1 ] == 112 ); //std::cout << ts_selector[ 1 ] << std::endl;
		TS_ASSERT( ts_selector[ 2 ] == 104 ); //std::cout << ts_selector[ 2 ] << std::endl;
		TS_ASSERT( ts_selector[ 3 ] == 150 ); //std::cout << ts_selector[ 3 ] << std::endl;
		TS_ASSERT( ts_selector[ 4 ] == 103 ); //std::cout << ts_selector[ 4 ] << std::endl;
		TS_ASSERT( ts_selector[ 5 ] == -101 ); //std::cout << ts_selector[ 5 ] << std::endl;

	}


};


