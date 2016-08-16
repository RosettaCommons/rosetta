// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/vector0.cxxtest.hh
/// @brief  vector0.test: test suite for utility::vector0
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/integer_mapping.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

#include <platform/types.hh>

#include <iostream>

using namespace utility;

class integer_mappingTests : public CxxTest::TestSuite {

public:
	typedef platform::Size          Size;
	typedef utility::subset_mapping subset_mapping;

public:

	/// @brief Subset_mapping should be constructed with 0 elements in the source
	/// enumeration and 0 elements in the destination enumeration
	void test_subset_mapping_default_ctor()
	{
		subset_mapping sm;
		TS_ASSERT( sm.source_size() == 0 );
		TS_ASSERT( sm.destination_size() == 0 );
	}

	void test_subset_mapping_size_ctor()
	{
		subset_mapping sm( 10 );
		TS_ASSERT( sm.source_size() == 10 );
		TS_ASSERT( sm.destination_size() == 0 );
		for ( Size ii = 1; ii <= 10; ++ii ) {
			TS_ASSERT( sm.s2d( ii ) == subset_mapping::UNMAPPED );
		}
	}

	void test_subset_mapping_odds_map() {
		subset_mapping sm( 10 );
		Size count = 0;
		for ( Size ii = 1; ii <= 10; ii += 2 ) {
			sm.set_next_correspondence( ii );
			++count;
			TS_ASSERT( sm.source_size() == 10 );
			TS_ASSERT( sm.destination_size() == count );
		}
		count = 0;
		for ( Size ii = 1; ii <= 10; ++ii ) {
			if ( ii % 2 == 0 ) {
				TS_ASSERT( sm.s2d( ii ) == subset_mapping::UNMAPPED );
			} else {
				++count;
				TS_ASSERT( sm.s2d( ii ) == count );
				TS_ASSERT( sm.d2s( count ) == ii );
			}
		}
	}

	void test_subset_mapping_copy_ctor()
	{
		subset_mapping sm( 10 );
		for ( Size ii = 1; ii <= 10; ii += 2 ) {
			sm.set_next_correspondence( ii );
		}
		subset_mapping sm2( sm );
		/// 1. Test that the mapping copied.
		Size count( 0 );
		TS_ASSERT( sm2.source_size() == 10 );
		TS_ASSERT( sm2.destination_size() == 5 );
		for ( Size ii = 1; ii <= 10; ++ii ) {
			if ( ii % 2 == 0 ) {
				TS_ASSERT( sm2.s2d( ii ) == subset_mapping::UNMAPPED );
			} else {
				++count;
				TS_ASSERT( sm2.s2d( ii ) == count );
				TS_ASSERT( sm2.d2s( count ) == ii );
			}
		}
		/// 2. Test that it was a deep copy -- modify the original, make sure the copy has remained the same.
		sm.set_next_correspondence( 2 );
		TS_ASSERT( sm.destination_size() == 6 );

		TS_ASSERT( sm2.source_size() == 10 );
		TS_ASSERT( sm2.destination_size() == 5 );

		count = 0;
		for ( Size ii = 1; ii <= 10; ++ii ) {
			if ( ii % 2 == 0 ) {
				TS_ASSERT( sm2.s2d( ii ) == subset_mapping::UNMAPPED );
			} else {
				++count;
				TS_ASSERT( sm2.s2d( ii ) == count );
				TS_ASSERT( sm2.d2s( count ) == ii );
			}
		}

	}

	void test_subset_mapping_out_of_bounds_set_next_coorespondence() {
		subset_mapping sm( 10 );
		try {
			sm.set_next_correspondence( 11 );
			TS_ASSERT( false ); /// should not reach here.
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_subset_mapping_out_of_bounds_set_next_coorespondence" << std::endl;
			//std::cout << e.msg() << std::endl;
			TS_ASSERT( e.msg() == "subset_mapping::set_next_correspondence recieved an out-of-bounds source id (11) with a source-enumeration size of 10" );
		}
	}

	void test_subset_mapping_overwrite_correspondence() {
		subset_mapping sm( 10 );
		sm.set_next_correspondence( 5 );
		try {
			sm.set_next_correspondence( 5 );
			TS_ASSERT( false ); /// should not reach here.
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cout << "test_subset_mapping_overwrite_correspondence" << std::endl;
			//std::cout << e.msg() << std::endl;
			TS_ASSERT( e.msg() == "subset_mapping::set_next_correspondence recieved an already-mapped source id (5) which had been previously assigned to destination id 1" );
		}
	}

}; // class heapTests

