// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/LexicographicalIterator.cxxtest.hh
/// @brief  test suite for utility::LexicographicalIterator
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Package headers
#include <utility/LexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>
#include <cxxtest/TestSuite.h>

/// C++ headers
#include <iostream>


class LexicographicalIteratorTests : public CxxTest::TestSuite {

private:
		utility::vector1< unsigned int > dimsizes_;
		utility::fixedsizearray1< platform::Size, 3 > fs_dimsizes_;

public:

	void setUp() {
		dimsizes_.resize( 3 );
		dimsizes_[ 1 ] = 4;
		dimsizes_[ 2 ] = 5;
		dimsizes_[ 3 ] = 3;
		std::copy( dimsizes_.begin(), dimsizes_.end(), fs_dimsizes_.begin() );
	}

	/// @brief Constructor Test
	void test_LexicographicalIterator_constructor() {

		utility::LexicographicalIterator li1( dimsizes_ );
		utility::LexicographicalIterator li2;

		TS_ASSERT( li1.dimsize( 1 ) == 4 );
		TS_ASSERT( li1.dimsize( 2 ) == 5 );
		TS_ASSERT( li1.dimsize( 3 ) == 3 );
	}

	void test_LexicographicalIterator_num_states_total() {

		utility::LexicographicalIterator li1( dimsizes_ );
		TS_ASSERT( li1.num_states_total() == 60 );

	}

	void test_LexicographicalIterator_increment() {
		utility::LexicographicalIterator li( dimsizes_ );
		li.begin();
		TS_ASSERT( li[ 1 ] == 1 && li[ 2 ] == 1 && li[ 3 ] == 1 );
		++li;
		TS_ASSERT( li[ 1 ] == 1 && li[ 2 ] == 1 && li[ 3 ] == 2 );
		++li; ++li;
		TS_ASSERT( li[ 1 ] == 1 && li[ 2 ] == 2 && li[ 3 ] == 1 );
	}

	void test_LexicographicalIterator_at_end() {
		utility::LexicographicalIterator li( dimsizes_ );
		li.begin();
		for ( unsigned int ii = 1; ii <= 60; ++ii ) {
			TS_ASSERT( ! li.at_end() );
			++li;
		}
		TS_ASSERT( li.at_end() );
	}

	void test_LexicographicalIterator_reset_state() {
		utility::LexicographicalIterator li( dimsizes_ );
		li.begin();
		for ( unsigned int ii = 1; ii <= 23; ++ii ) ++li;
		unsigned int d1( li[ 1 ]), d2( li[ 2 ] ), d3( li[ 3 ] );
		///std::cout << "d1 " << d1 << " d2 " << d2 << " d3 " << d3 << std::endl;

		unsigned int state = li.index();
		///std::cout << "state: " << state << std::endl;
		++li; ++li; ++li; ++li; ++li; ++li;
		li.set_position_from_index( state );

		///std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;

		TS_ASSERT( d1 == li[ 1 ] );
		TS_ASSERT( d2 == li[ 2 ] );
		TS_ASSERT( d3 == li[ 3 ] );
		TS_ASSERT( ! li.at_end() );

		/// Test that the last state sets all the dimensions at their maximal values.
		li.set_position_from_index( 60 );
		TS_ASSERT( li[ 1 ] == 4 );
		TS_ASSERT( li[ 2 ] == 5 );
		TS_ASSERT( li[ 3 ] == 3 );
		TS_ASSERT( !li.at_end() );

		/// Test that initializing from an index past the boundary triggers at_end().
		li.set_position_from_index( 80 );
		TS_ASSERT( li.at_end() );
	}

	void test_LexicographicalIterator_continue_at_dim() {
		utility::LexicographicalIterator li( dimsizes_ );
		li.begin();
		for ( unsigned int ii = 1; ii <= 5; ++ii ) ++li;

		TS_ASSERT( li[ 1 ] == 1 );
		TS_ASSERT( li[ 2 ] == 2 );
		TS_ASSERT( li[ 3 ] == 3 );

		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;

		unsigned int const state( li.index() );

		li.continue_at_dimension( 3 );
		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;
		TS_ASSERT( li[ 1 ] == 1 );
		TS_ASSERT( li[ 2 ] == 3 );
		TS_ASSERT( li[ 3 ] == 1 );

		li.set_position_from_index( state ); // restore

		li.continue_at_dimension( 2 );
		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;
		TS_ASSERT( li[ 1 ] == 1 );
		TS_ASSERT( li[ 2 ] == 3 );
		TS_ASSERT( li[ 3 ] == 1 );

		li.set_position_from_index( state ); // restore

		li.continue_at_dimension( 1 );
		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;
		TS_ASSERT( li[ 1 ] == 2 );
		TS_ASSERT( li[ 2 ] == 1 );
		TS_ASSERT( li[ 3 ] == 1 );
	}

//// FIXED SIZE LEX

	/// @brief Constructor Test
	void test_FixedSizeLexicographicalIterator_constructor() {

		utility::FixedSizeLexicographicalIterator< 3 > li1( fs_dimsizes_ );

		TS_ASSERT( li1.dimsize( 1 ) == 4 );
		TS_ASSERT( li1.dimsize( 2 ) == 5 );
		TS_ASSERT( li1.dimsize( 3 ) == 3 );
	}

	/// @brief Constructor Test
	void test_FixedSizeLexicographicalIterator_set_dim_sizes() {

		utility::FixedSizeLexicographicalIterator< 3 > li1;
		li1.set_dimension_sizes( fs_dimsizes_ );

		TS_ASSERT( li1.dimsize( 1 ) == 4 );
		TS_ASSERT( li1.dimsize( 2 ) == 5 );
		TS_ASSERT( li1.dimsize( 3 ) == 3 );

		TS_ASSERT( li1[ 1 ] == 1 );
		TS_ASSERT( li1[ 2 ] == 1 );
		TS_ASSERT( li1[ 3 ] == 1 );
	}

	void test_FixedSizeLexicographicalIterator_num_states_total() {

		utility::FixedSizeLexicographicalIterator< 3 > li1( fs_dimsizes_ );
		TS_ASSERT( li1.num_states_total() == 60 );

	}

	void test_FixedSizeLexicographicalIterator_increment() {
		utility::FixedSizeLexicographicalIterator< 3 > li( fs_dimsizes_ );
		li.begin();
		TS_ASSERT( li[ 1 ] == 1 && li[ 2 ] == 1 && li[ 3 ] == 1 );
		++li;
		TS_ASSERT( li[ 1 ] == 1 && li[ 2 ] == 1 && li[ 3 ] == 2 );
		++li; ++li;
		TS_ASSERT( li[ 1 ] == 1 && li[ 2 ] == 2 && li[ 3 ] == 1 );
	}

	void test_FixedSizeLexicographicalIterator_at_end() {
		utility::FixedSizeLexicographicalIterator< 3 > li( fs_dimsizes_ );
		li.begin();
		for ( unsigned int ii = 1; ii <= 60; ++ii ) {
			TS_ASSERT( ! li.at_end() );
			++li;
		}
		TS_ASSERT( li.at_end() );
	}

	void test_FixedSizeLexicographicalIterator_reset_state() {
		utility::FixedSizeLexicographicalIterator< 3 > li( fs_dimsizes_ );
		li.begin();
		for ( unsigned int ii = 1; ii <= 23; ++ii ) ++li;
		unsigned int d1( li[ 1 ]), d2( li[ 2 ] ), d3( li[ 3 ] );
		///std::cout << "d1 " << d1 << " d2 " << d2 << " d3 " << d3 << std::endl;

		unsigned int state = li.index();
		///std::cout << "state: " << state << std::endl;
		++li; ++li; ++li; ++li; ++li; ++li;
		li.set_position_from_index( state );

		///std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;

		TS_ASSERT( d1 == li[ 1 ] );
		TS_ASSERT( d2 == li[ 2 ] );
		TS_ASSERT( d3 == li[ 3 ] );
		TS_ASSERT( ! li.at_end() );

		/// Test that the last state sets all the dimensions at their maximal values.
		li.set_position_from_index( 60 );
		TS_ASSERT( li[ 1 ] == 4 );
		TS_ASSERT( li[ 2 ] == 5 );
		TS_ASSERT( li[ 3 ] == 3 );
		TS_ASSERT( !li.at_end() );

		/// Test that initializing from an index past the boundary triggers at_end().
		li.set_position_from_index( 80 );
		TS_ASSERT( li.at_end() );
	}

	void test_FixedSizeLexicographicalIterator_continue_at_dim() {
		utility::FixedSizeLexicographicalIterator< 3 > li( fs_dimsizes_ );
		li.begin();
		for ( unsigned int ii = 1; ii <= 5; ++ii ) ++li;

		TS_ASSERT( li[ 1 ] == 1 );
		TS_ASSERT( li[ 2 ] == 2 );
		TS_ASSERT( li[ 3 ] == 3 );

		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;

		unsigned int const state( li.index() );

		li.continue_at_dimension( 3 );
		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;
		TS_ASSERT( li[ 1 ] == 1 );
		TS_ASSERT( li[ 2 ] == 3 );
		TS_ASSERT( li[ 3 ] == 1 );

		li.set_position_from_index( state ); // restore

		li.continue_at_dimension( 2 );
		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;
		TS_ASSERT( li[ 1 ] == 1 );
		TS_ASSERT( li[ 2 ] == 3 );
		TS_ASSERT( li[ 3 ] == 1 );

		li.set_position_from_index( state ); // restore

		li.continue_at_dimension( 1 );
		//std::cout << "li[ 1 ] " << li[ 1 ] << " li[ 2 ] " << li[ 2 ] << " li[ 3 ] " << li[ 3 ] << std::endl;
		TS_ASSERT( li[ 1 ] == 2 );
		TS_ASSERT( li[ 2 ] == 1 );
		TS_ASSERT( li[ 3 ] == 1 );
	}


};


