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
#include <utility/recent_history_queue.hh>

#include <iostream>

using namespace utility;

class RecentHistoryQueueTests : public CxxTest::TestSuite {

public:

	/// @brief very basic -- initialize a list in a particular order and make sure
	/// that the list reflects that order
	void test_recent_history_queue_fill_queue()
	{
		recent_history_queue rhq( 10, 5 );
		TS_ASSERT( rhq.num_elements() ==  10 );
		TS_ASSERT( rhq.history_size() == 5 );
		TS_ASSERT( rhq.curr_num_in_recent_history() == 0 );

		int bumped = -1;
		bumped = rhq.push_to_front_of_history_queue( 4 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.curr_num_in_recent_history() == 1 );
		TS_ASSERT( rhq.pos_in_history_queue( 4 ) == 1 );
		TS_ASSERT( rhq.pos_in_history_queue( 3 ) == 0 );

		bumped = rhq.push_to_front_of_history_queue( 3 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.curr_num_in_recent_history() == 2 );
		TS_ASSERT( rhq.pos_in_history_queue( 3 ) == 2 );
		TS_ASSERT( rhq.pos_in_history_queue( 7 ) == 0 );

		bumped = rhq.push_to_front_of_history_queue( 7 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.curr_num_in_recent_history() == 3 );
		TS_ASSERT( rhq.pos_in_history_queue( 7 ) == 3 );
		TS_ASSERT( rhq.pos_in_history_queue( 1 ) == 0 );

		bumped = rhq.push_to_front_of_history_queue( 1 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.curr_num_in_recent_history() == 4 );
		TS_ASSERT( rhq.pos_in_history_queue( 1 ) == 4 );
		TS_ASSERT( rhq.pos_in_history_queue( 2 ) == 0 );

		bumped = rhq.push_to_front_of_history_queue( 2 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.curr_num_in_recent_history() == 5 );
		TS_ASSERT( rhq.pos_in_history_queue( 2 ) == 5 );

	}

	void test_reset_recent_history_queue_tail() {
		recent_history_queue rhq( 10, 5 );
		rhq.push_to_front_of_history_queue( 4 );
		rhq.push_to_front_of_history_queue( 3 );
		rhq.push_to_front_of_history_queue( 7 );
		rhq.push_to_front_of_history_queue( 1 );
		rhq.push_to_front_of_history_queue( 2 );

		TS_ASSERT( rhq.head_of_queue() == 5 );
		TS_ASSERT( rhq.end_of_queue() == 1 );

		int bumped = rhq.push_to_front_of_history_queue( 4 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.head_of_queue() == 1 );
		TS_ASSERT( rhq.end_of_queue() == 2 );
	}

	void test_reset_recent_history_queue_head() {
		recent_history_queue rhq( 10, 5 );
		rhq.push_to_front_of_history_queue( 4 );
		rhq.push_to_front_of_history_queue( 3 );
		rhq.push_to_front_of_history_queue( 7 );
		rhq.push_to_front_of_history_queue( 1 );
		rhq.push_to_front_of_history_queue( 2 );

		TS_ASSERT( rhq.head_of_queue() == 5 );
		TS_ASSERT( rhq.end_of_queue() == 1 );

		int bumped = rhq.push_to_front_of_history_queue( 2 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.head_of_queue() == 5 );
		TS_ASSERT( rhq.end_of_queue() == 1 );
	}

	void test_reset_recent_history_queue_middle() {
		recent_history_queue rhq( 10, 5 );
		rhq.push_to_front_of_history_queue( 4 );
		rhq.push_to_front_of_history_queue( 3 );
		rhq.push_to_front_of_history_queue( 7 );
		rhq.push_to_front_of_history_queue( 1 );
		rhq.push_to_front_of_history_queue( 2 );

		TS_ASSERT( rhq.head_of_queue() == 5 );
		TS_ASSERT( rhq.end_of_queue() == 1 );

		int bumped = rhq.push_to_front_of_history_queue( 7 );
		TS_ASSERT( bumped == 0 );
		TS_ASSERT( rhq.head_of_queue() == 3 );
		TS_ASSERT( rhq.end_of_queue() == 1 );
	}

	void test_bump_an_element_from_the_history() {
		recent_history_queue rhq( 10, 5 );
		rhq.push_to_front_of_history_queue( 4 );
		rhq.push_to_front_of_history_queue( 3 );
		rhq.push_to_front_of_history_queue( 7 );
		rhq.push_to_front_of_history_queue( 1 );
		rhq.push_to_front_of_history_queue( 2 );

		int bumped = rhq.push_to_front_of_history_queue( 5 );
		TS_ASSERT( bumped == 1 );
		TS_ASSERT( rhq.pos_in_history_queue( 4 ) == 0 );
		TS_ASSERT( rhq.pos_in_history_queue( 5 ) == 1 );
	}

	void test_bump_an_element_from_the_history_after_resent_queue_end() {
		recent_history_queue rhq( 10, 5 );
		rhq.push_to_front_of_history_queue( 4 );
		rhq.push_to_front_of_history_queue( 3 );
		rhq.push_to_front_of_history_queue( 7 );
		rhq.push_to_front_of_history_queue( 1 );
		rhq.push_to_front_of_history_queue( 2 );

		rhq.push_to_front_of_history_queue( 4 ); // push 4 to the front

		int bumped = rhq.push_to_front_of_history_queue( 5 );
		TS_ASSERT( bumped == 2 );
		TS_ASSERT( rhq.pos_in_history_queue( 4 ) == 1 );
		TS_ASSERT( rhq.pos_in_history_queue( 3 ) == 0 );
		TS_ASSERT( rhq.pos_in_history_queue( 5 ) == 2 );
	}

}; // class RecentHistoryQueueTests

