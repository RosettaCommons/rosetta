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
#include <utility/OrderedTuple.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>

#include <cxxtest/TestSuite.h>

/// C++ headers
#include <iostream>
#include <list>

using namespace utility;

class OrderedTupleTests : public CxxTest::TestSuite {

public:
	typedef fixedsizearray1< int, 3 > int3;
	typedef OrderedTuple< int3 > int3tuple;

public:

	void setUp() {
	}

	/// @brief Test that the OrderedTuple works for a fixedsizearray1
	void test_OrderedTuple_data_assignment_fsa() {
		int3 array; array[ 1 ] = 1; array[ 2 ] = 15; array[ 3 ] = 12;
		int3tuple tuple;
		tuple.assign_data( array );

		int3tuple::const_iterator iter( tuple.begin() ), iter_end( tuple.end() );
		TS_ASSERT( *iter == 1 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( *iter == 15 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( *iter == 12 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( iter == iter_end );
	}

	/// @brief Test that the ordered tuple works for a utility::vector1
	void test_OrderedTuple_data_assignment_v1() {
		utility::vector1< int > array( 3 );
		array[ 1 ] = 1; array[ 2 ] = 15; array[ 3 ] = 12;
		OrderedTuple< utility::vector1< int > > tuple;
		tuple.assign_data( array );

		OrderedTuple< utility::vector1< int > >::const_iterator
			iter( tuple.begin() ), iter_end( tuple.end() );
		TS_ASSERT( *iter == 1 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( *iter == 15 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( *iter == 12 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( iter == iter_end );
	}

	/// @brief Test that the ordered tuple works for a list
	void test_OrderedTuple_data_assignment_list() {
		std::list< int > list3;
		list3.push_back(  1 ); list3.push_back( 15 ); list3.push_back( 12 );
		OrderedTuple< std::list< int > > tuple;
		tuple.assign_data( list3 );

		OrderedTuple< OrderedTuple< std::list< int > > >::const_iterator
			iter( tuple.begin() ), iter_end( tuple.end() );
		TS_ASSERT( *iter == 1 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( *iter == 15 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( *iter == 12 );
		TS_ASSERT( iter != iter_end );
		++iter;
		TS_ASSERT( iter == iter_end );
	}

	void test_OrderedTuple_compare_fsa() {
		int3 arrayA; arrayA[ 1 ] = 1; arrayA[ 2 ] = 15; arrayA[ 3 ] = 12;
		int3tuple tupleA;
		tupleA.assign_data( arrayA );

		int3 arrayB; arrayB[ 1 ] = 1; arrayB[ 2 ] = 15; arrayB[ 3 ] = 11;
		int3tuple tupleB;
		tupleB.assign_data( arrayB );

		int3 arrayC; arrayC[ 1 ] = 1; arrayC[ 2 ] = 3; arrayC[ 3 ] = 12;
		int3tuple tupleC;
		tupleC.assign_data( arrayC );

		TS_ASSERT( tupleB < tupleA );
		TS_ASSERT( tupleC < tupleB );
		TS_ASSERT( tupleC < tupleA );

	}

};


