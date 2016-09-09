// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/ArrayPool.cxxtest.hh
/// @brief  test suite for utility::graph::ArrayPool.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <utility/graph/ArrayPool.hh>

//Auto Headers
#include <utility/vector1.hh>

//Types
#include <platform/types.hh>


using namespace platform;
using namespace utility;
using namespace utility::graph;


// --------------- Test Class --------------- //

class ArrayPoolTests : public CxxTest::TestSuite {

public:


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
		// Being a smart pointer, g should be destructed and "free'd" correctly, but a g->delete_everything()
		// could be placed here, if desired.
	}


	// --------------- Test Cases --------------- //
	void test_pool_default_ctor() {
		ArrayPool< double > pool;
		TS_ASSERT( pool.empty() );
		TS_ASSERT( pool.array_size() == 0 );
		TS_ASSERT( pool.block_size() == 32 );
	}

	void test_pool_1argument_ctor() {
		ArrayPool< double > pool( 256 );
		TS_ASSERT( pool.empty() );
		TS_ASSERT( pool.array_size() == 0 );
		TS_ASSERT( pool.block_size() == 256 );
	}

	void test_pool_2argument_ctor() {
		ArrayPool< double > pool( 256, 5 );
		TS_ASSERT( pool.empty() );
		TS_ASSERT( pool.array_size() == 5 );
		TS_ASSERT( pool.block_size() == 256 );
	}

	void test_allocate_one_array() {
		ArrayPool< double > pool( 32, 5 );
		ArrayPoolElement< double > element = pool.new_array();
		TS_ASSERT( pool.empty() == false );

		element[ 0 ] = 5;
		element[ 1 ] = 12;
		element[ 2 ] = 0.5;
		element[ 3 ] = 0.25;
		element[ 4 ] = 100;
		TS_ASSERT( element[ 0 ] == 5 );
		TS_ASSERT( element[ 1 ] == 12 );
		TS_ASSERT( element[ 2 ] == 0.5 );
		TS_ASSERT( element[ 3 ] == 0.25 );
		TS_ASSERT( element[ 4 ] == 100 );


		pool.deallocate_array( element );
		TS_ASSERT( pool.empty() );
	}

	void test_allocate_one_block() {
		ArrayPool< double > pool( 8, 5 );
		std::vector< ArrayPoolElement< double > > arrays;
		arrays.reserve( 8 );
		for ( platform::Size ii = 0; ii < 8; ++ii ) {
			arrays.push_back( pool.new_array() );
			TS_ASSERT( arrays[ ii ].valid() );
		}
		TS_ASSERT( pool.nblocks() == 1 );
		TS_ASSERT( pool.noutstanding() == 8 );
		for ( platform::Size ii = 0; ii < 8; ++ii ) {
			TS_ASSERT( arrays[ ii ].valid() );
			pool.deallocate_array( arrays[ ii ] );
			TS_ASSERT( ! arrays[ ii ].valid() );
		}
	}

	void test_allocate_two_blocks() {
		ArrayPool< double > pool( 8, 5 );
		std::vector< ArrayPoolElement< double > > arrays;
		arrays.reserve( 9 );
		for ( platform::Size ii = 0; ii < 8; ++ii ) {
			arrays.push_back( pool.new_array() );
			TS_ASSERT( arrays[ ii ].valid() );
		}
		TS_ASSERT( pool.nblocks() == 1 );
		TS_ASSERT( pool.noutstanding() == 8 );

		arrays.push_back( pool.new_array() );

		TS_ASSERT( pool.nblocks() == 2 );
		TS_ASSERT( pool.noutstanding() == 9 );
		for ( platform::Size ii = 0; ii < 9; ++ii ) {
			TS_ASSERT( arrays[ ii ].valid() );
			pool.deallocate_array( arrays[ ii ] );
			TS_ASSERT( ! arrays[ ii ].valid() );
		}
	}

	void test_reuse_space() {
		ArrayPool< double > pool( 8, 5 );
		std::vector< ArrayPoolElement< double > > arrays;
		arrays.reserve( 8 );
		for ( platform::Size ii = 0; ii < 8; ++ii ) {
			arrays.push_back( pool.new_array() );
		}
		TS_ASSERT( pool.noutstanding() == 8 );
		TS_ASSERT( pool.nblocks() == 1 );

		pool.deallocate_array( arrays[ 5 ] );
		TS_ASSERT( pool.noutstanding() == 7 );

		arrays[ 5 ] = pool.new_array();
		TS_ASSERT( pool.noutstanding() == 8 );
		TS_ASSERT( pool.nblocks() == 1 );

		arrays.push_back( pool.new_array() );
		TS_ASSERT( pool.noutstanding() == 9 );
		TS_ASSERT( pool.nblocks() == 2 );

		for ( platform::Size ii = 0; ii < 9; ++ii ) {
			TS_ASSERT( arrays[ ii ].valid() );
			pool.deallocate_array( arrays[ ii ] );
			TS_ASSERT( ! arrays[ ii ].valid() );
		}

	}


};


