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
#include <utility/heap.hh>

#include <iostream>

using namespace utility;

class heapTests : public CxxTest::TestSuite {

public:

	/// @brief Size + Value Constructor Test
	void test_heap_constructor() {

		heap h( 5 );
		TS_ASSERT( h.size() == 0 );
		TS_ASSERT( h.capacity() == 5 );
	}

	void test_heap_insert_one_element() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1, 0.5, err );
		TS_ASSERT( ! err );
		TS_ASSERT( h.size() == 1 );
		TS_ASSERT( h.capacity() == 5 );
	}

	void test_heap_remove_element_from_one_element_heap() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1, 0.5, err );

		int val; float coval;
		h.heap_extract( val, coval, err );
		TS_ASSERT( ! err );
		TS_ASSERT( val == 1 );
		TS_ASSERT( coval == 0.5 );
		TS_ASSERT( h.size() == 0 );
		TS_ASSERT( h.capacity() == 5 );
	}

	void test_heap_add_two_elememnts() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1, 0.5, err );
		TS_ASSERT( ! err );
		h.heap_insert( 2, 0.25, err );
		TS_ASSERT( ! err );
		TS_ASSERT( h.size() == 2 );
		TS_ASSERT( h.capacity() == 5 );
	}


	void test_heap_remove_element_from_two_element_heap() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1, 0.5, err );
		h.heap_insert( 2, 0.25, err );

		int val; float coval;
		h.heap_extract( val, coval, err );
		TS_ASSERT( ! err );
		TS_ASSERT( val == 2 );
		TS_ASSERT( coval == 0.25 );
	}

	void test_heap_exceed_capacity() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1,  0.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 2,  0.25, err ); TS_ASSERT( ! err );
		h.heap_insert( 3,  0.75, err ); TS_ASSERT( ! err );
		h.heap_insert( 4,  1.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 5, -0.5,  err ); TS_ASSERT( ! err );

		h.heap_insert( 6, 3.5, err ); TS_ASSERT( err );

		TS_ASSERT( h.coval(0) == 0.25 );
		TS_ASSERT( h.coval(1) == 0.5 );
		TS_ASSERT( h.coval(2) == 0.75 );
		TS_ASSERT( h.coval(3) == 1.5 );
		TS_ASSERT( h.coval(4) == 3.5 );

	}

	void test_heap_remove_from_full_heap() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1,  0.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 2,  0.25, err ); TS_ASSERT( ! err );
		h.heap_insert( 3,  0.75, err ); TS_ASSERT( ! err );
		h.heap_insert( 4,  1.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 5, -0.5,  err ); TS_ASSERT( ! err );

		int val; float coval;
		h.heap_extract( val, coval, err );
		TS_ASSERT( ! err );
		TS_ASSERT( val == 5 );
		TS_ASSERT( coval == -0.5 );
		TS_ASSERT( h.size() == 4 );

		TS_ASSERT( h.coval(0) == 0.25 );
		TS_ASSERT( h.coval(1) == 0.5 );
		TS_ASSERT( h.coval(2) == 0.75 );
		TS_ASSERT( h.coval(3) == 1.5 );
	}


	void test_heap_decrease_coval() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1,  0.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 2,  0.25, err ); TS_ASSERT( ! err );
		h.heap_insert( 3,  0.75, err ); TS_ASSERT( ! err );
		h.heap_insert( 4,  1.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 5, -0.5,  err ); TS_ASSERT( ! err );

		h.reset_coval( 3, -0.25 );
		//for ( int ii = 0; ii < 5; ++ii ) std::cout << "TS_ASSERT( h.coval(" << ii << ") == " << h.coval( ii ) << " );" << std::endl;
		//for ( int ii = 0; ii < 5; ++ii ) std::cout << "TS_ASSERT( h.val(" << ii << ") == " << h.val( ii ) << " );" << std::endl;
		TS_ASSERT( h.coval(0) == -0.5 );
		TS_ASSERT( h.coval(1) == 0.25 );
		TS_ASSERT( h.coval(2) == -0.25 );
		TS_ASSERT( h.coval(3) == 1.5 );
		TS_ASSERT( h.coval(4) == 0.5 );
		TS_ASSERT( h.val(0) == 5 );
		TS_ASSERT( h.val(1) == 2 );
		TS_ASSERT( h.val(2) == 3 );
		TS_ASSERT( h.val(3) == 4 );
		TS_ASSERT( h.val(4) == 1 );

		h.reset_coval( 3, -0.75 );
		//for ( int ii = 0; ii < 5; ++ii ) std::cout << "TS_ASSERT( h.coval(" << ii << ") == " << h.coval( ii ) << " );" << std::endl;
		//for ( int ii = 0; ii < 5; ++ii ) std::cout << "TS_ASSERT( h.val(" << ii << ") == " << h.val( ii ) << " );" << std::endl;
		TS_ASSERT( h.coval(0) == -0.75 );
		TS_ASSERT( h.coval(1) == 0.25 );
		TS_ASSERT( h.coval(2) == -0.5 );
		TS_ASSERT( h.coval(3) == 1.5 );
		TS_ASSERT( h.coval(4) == 0.5 );
		TS_ASSERT( h.val(0) == 3 );
		TS_ASSERT( h.val(1) == 2 );
		TS_ASSERT( h.val(2) == 5 );
		TS_ASSERT( h.val(3) == 4 );
		TS_ASSERT( h.val(4) == 1 );

	}

	void test_heap_increase_coval() {
		heap h( 5 );
		bool err;
		h.heap_insert( 1,  0.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 2,  0.25, err ); TS_ASSERT( ! err );
		h.heap_insert( 3,  0.75, err ); TS_ASSERT( ! err );
		h.heap_insert( 4,  1.5,  err ); TS_ASSERT( ! err );
		h.heap_insert( 5, -0.5,  err ); TS_ASSERT( ! err );

		h.reset_coval( 1, 1.25 );
		//for ( int ii = 0; ii < 5; ++ii ) std::cout << "TS_ASSERT( h.coval(" << ii << ") == " << h.coval( ii ) << " );" << std::endl;
		//for ( int ii = 0; ii < 5; ++ii ) std::cout << "TS_ASSERT( h.val(" << ii << ") == " << h.val( ii ) << " );" << std::endl;
		TS_ASSERT( h.coval(0) == -0.5 );
		TS_ASSERT( h.coval(1) == 0.25 );
		TS_ASSERT( h.coval(2) == 0.75 );
		TS_ASSERT( h.coval(3) == 1.5 );
		TS_ASSERT( h.coval(4) == 1.25 );
		TS_ASSERT( h.val(0) == 5 );
		TS_ASSERT( h.val(1) == 2 );
		TS_ASSERT( h.val(2) == 3 );
		TS_ASSERT( h.val(3) == 4 );
		TS_ASSERT( h.val(4) == 1 );
	}

}; // class heapTests

