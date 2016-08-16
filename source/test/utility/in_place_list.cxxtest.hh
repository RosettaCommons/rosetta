// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/in_place_list.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/in_place_list.hh>

#include <iostream>

using namespace utility;

class in_place_listTests : public CxxTest::TestSuite {

public:

	/// @brief very basic -- initialize a list in a particular order and make sure
	/// that the list reflects that order
	void test_create_in_place_list() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 1 );
					break;
				case 2 :
					TS_ASSERT( ind == 3 );
					break;
				case 3 :
					TS_ASSERT( ind == 5 );
					break;
				case 4 :
					TS_ASSERT( ind == 2 );
					break;
				case 5 :
					TS_ASSERT( ind == 4 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}
		TS_ASSERT( count == 6 );
		TS_ASSERT( intlist.head() == 1 );
		TS_ASSERT( intlist.tail() == 4 );
	}

	/// @brief test the ability to drop a middle element from the list
	void test_drop_element_in_place_list() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		intlist.remove( 5 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 1 );
					break;
				case 2 :
					TS_ASSERT( ind == 3 );
					break;
				case 3 :
					TS_ASSERT( ind == 2 );
					break;
				case 4 :
					TS_ASSERT( ind == 4 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}
		TS_ASSERT( count == 5 );
		TS_ASSERT( ! intlist[ 5 ].in_list() );

	}

	/// @brief make sure you can remove the first element in a list
	void test_in_place_list_remove_from_head() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		intlist.remove( 1 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 3 );
					break;
				case 2 :
					TS_ASSERT( ind == 5 );
					break;
				case 3 :
					TS_ASSERT( ind == 2 );
					break;
				case 4 :
					TS_ASSERT( ind == 4 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}

		TS_ASSERT( count == 5 );
		TS_ASSERT( ! intlist[ 1 ].in_list() );
		TS_ASSERT( intlist.head() == 3 );
		TS_ASSERT( intlist.tail() == 4 );

	}

	/// @brief make sure you can remove the last element in the list
	void test_in_place_list_remove_from_tail() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		intlist.remove( 4 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 1 );
					break;
				case 2 :
					TS_ASSERT( ind == 3 );
					break;
				case 3 :
					TS_ASSERT( ind == 5 );
					break;
				case 4 :
					TS_ASSERT( ind == 2 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}

		TS_ASSERT( count == 5 );
		TS_ASSERT( ! intlist[ 4 ].in_list() );
		TS_ASSERT( intlist.head() == 1 );
		TS_ASSERT( intlist.tail() == 2 );

	}

	void test_in_place_list_move_to_front_from_middle() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );


		intlist.move_to_front( 5 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 5 );
					break;
				case 2 :
					TS_ASSERT( ind == 1 );
					break;
				case 3 :
					TS_ASSERT( ind == 3 );
					break;
				case 4 :
					TS_ASSERT( ind == 2 );
					break;
				case 5 :
					TS_ASSERT( ind == 4 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}

		TS_ASSERT( count == 6 );
		TS_ASSERT( intlist.head() == 5 );
		TS_ASSERT( intlist.tail() == 4 );

	}

	void test_in_place_list_move_to_front_from_back() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );


		intlist.move_to_front( 4 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 4 );
					break;
				case 2 :
					TS_ASSERT( ind == 1 );
					break;
				case 3 :
					TS_ASSERT( ind == 3 );
					break;
				case 4 :
					TS_ASSERT( ind == 5 );
					break;
				case 5 :
					TS_ASSERT( ind == 2 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}

		TS_ASSERT( count == 6 );
		TS_ASSERT( intlist.head() == 4 );
		TS_ASSERT( intlist.tail() == 2 );

	}

	void test_in_place_list_move_to_front_from_front() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		intlist.move_to_front( 1 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 1 );
					break;
				case 2 :
					TS_ASSERT( ind == 3 );
					break;
				case 3 :
					TS_ASSERT( ind == 5 );
					break;
				case 4 :
					TS_ASSERT( ind == 2 );
					break;
				case 5 :
					TS_ASSERT( ind == 4 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}
		TS_ASSERT( count == 6 );
		TS_ASSERT( intlist.head() == 1 );
		TS_ASSERT( intlist.tail() == 4 );
	}

	void test_in_place_list_remove_all_but_one() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		intlist.remove( 1 );
		intlist.remove( 3 );
		intlist.remove( 2 );
		intlist.remove( 4 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 5 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}
		TS_ASSERT( count == 2 );
		TS_ASSERT( ! intlist[ 1 ].in_list() );
		TS_ASSERT( ! intlist[ 2 ].in_list() );
		TS_ASSERT( ! intlist[ 3 ].in_list() );
		TS_ASSERT( ! intlist[ 4 ].in_list() );
		TS_ASSERT( intlist.head() == 5 );
		TS_ASSERT( intlist.tail() == 5 );
	}

	void test_in_place_list_remove_all_but_head() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		intlist.remove( 3 );
		intlist.remove( 2 );
		intlist.remove( 5 );
		intlist.remove( 4 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 1 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}
		TS_ASSERT( count == 2 );
		TS_ASSERT( ! intlist[ 2 ].in_list() );
		TS_ASSERT( ! intlist[ 3 ].in_list() );
		TS_ASSERT( ! intlist[ 4 ].in_list() );
		TS_ASSERT( ! intlist[ 5 ].in_list() );
		TS_ASSERT( intlist.head() == 1 );
		TS_ASSERT( intlist.tail() == 1 );
	}

	void test_in_place_list_remove_all_but_tail() {
		in_place_list< platform::Size > intlist( 5, 0 );

		for ( platform::Size ii = 1; ii <= intlist.size(); ++ii ) {
			intlist[ ii ].data() = ii;
		}

		intlist.move_to_front( 4 );
		intlist.move_to_front( 2 );
		intlist.move_to_front( 5 );
		intlist.move_to_front( 3 );
		intlist.move_to_front( 1 );

		intlist.remove( 3 );
		intlist.remove( 2 );
		intlist.remove( 5 );
		intlist.remove( 1 );

		platform::Size count = 1;
		for ( platform::Size ind = intlist.head(); ind != intlist.end(); ind = intlist[ ind ].next() ) {
			switch ( count ) {
				case 1 :
					TS_ASSERT( ind == 4 );
					break;
				default:
					TS_ASSERT( false );
			}
			TS_ASSERT( intlist[ ind ].in_list() );
			count += 1;
		}
		TS_ASSERT( count == 2 );
		TS_ASSERT( ! intlist[ 1 ].in_list() );
		TS_ASSERT( ! intlist[ 2 ].in_list() );
		TS_ASSERT( ! intlist[ 3 ].in_list() );
		TS_ASSERT( ! intlist[ 5 ].in_list() );
		TS_ASSERT( intlist.head() == 4 );
		TS_ASSERT( intlist.tail() == 4 );
	}


}; // class heapTests

