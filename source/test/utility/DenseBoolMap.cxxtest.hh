// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/DenseBoolMap.cxxtest.hh
/// @brief  test suite for utility::DenseBoolMap
/// @author Jack Maguire, jack@med.unc.edu

// Package headers
#include <utility/DenseBoolMap.hh>
#include <cxxtest/TestSuite.h>
#include <utility/vector1.hh>

class DenseBoolMapTests : public CxxTest::TestSuite {

public:
	void test_dense_bool_map(){

		/*
		Alright this is kind of silly but clang does not like
		TS_ASSERT_EQUALS( bool i, bool j ) when the compiler
		can determine j because the macro converts each argument
		to a string and clang can not convert "true" to the
		correct format of string. All of the values in the
		utility::vector1< bool >'s are known at compile time so
		we need to use TS_ASSERT( bool i ^ ! bool j ) to test
		for equality. I think TS_ASSERT( bool i == bool j ) is
		equivalent but someone on stack overflow was pessimistic
		that the == operator would always perform the way we want
		it to. Instead of looking up a second opinion or do any
		tests, I just used "^ !".

		-Jack Maguire
		*/

		utility::vector1< bool > target_settings1( 100, false );
		utility::vector1< bool > target_settings2( 100, false );
		utility::DenseBoolMap< 100, 1 > bool_map;

		for ( int i=1; i<101; ++i ) {
			target_settings1[ i ] = ( i % 5 == 0 );
			target_settings2[ i ] = ( i % 3 == 0 );
			bool_map.set( i, target_settings1[ i ] );
		}

		for ( int i=1; i<101; ++i ) {
			bool const bool_map_result = bool_map.get( i );
			TS_ASSERT( bool_map_result == target_settings1[ i ] );
		}

		TS_ASSERT( bool_map.get< 1 >() ^ ! target_settings1[ 1 ] );
		TS_ASSERT( bool_map.get< 2 >() ^ ! target_settings1[ 2 ] );
		TS_ASSERT( bool_map.get< 3 >() ^ ! target_settings1[ 3 ] );
		TS_ASSERT( bool_map.get< 4 >() ^ ! target_settings1[ 4 ] );
		TS_ASSERT( bool_map.get< 100 >() ^ ! target_settings1[ 100 ] );

		for ( int i=1; i<101; ++i ) {
			bool_map.set( i, target_settings2[ i ] );
		}

		for ( int i=1; i<101; ++i ) {
			bool const bool_map_result = bool_map.get( i );
			TS_ASSERT( bool_map_result == target_settings2[ i ] );
		}

		TS_ASSERT( bool_map.get< 1 >() ^ ! target_settings2[ 1 ] );
		TS_ASSERT( bool_map.get< 2 >() ^ ! target_settings2[ 2 ] );
		TS_ASSERT( bool_map.get< 3 >() ^ ! target_settings2[ 3 ] );
		TS_ASSERT( bool_map.get< 4 >() ^ ! target_settings2[ 4 ] );
		TS_ASSERT( bool_map.get< 100 >() ^ ! target_settings2[ 100 ] );

	}

};
