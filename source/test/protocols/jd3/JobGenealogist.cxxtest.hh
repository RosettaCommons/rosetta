// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/JobGenealogist.cxxtest.hh
/// @brief  test suite for the JobGenealogist
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/jd3/JobGenealogist.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>

using namespace protocols::jd3;

class JobGenealogistTests : public CxxTest::TestSuite
{
public:

	void setUp() {
	}

	void test_genealogist_set_target_num_jobs() {
		JobGenealogist jg;
		jg.set_num_nodes( 5 );
		jg.set_target_num_jobs_for_node( 1, 5 );
		jg.set_target_num_jobs_for_node( 2, 10 );
		jg.set_target_num_jobs_for_node( 4, 10 );
		jg.set_target_num_jobs_for_node( 3, 10 );
		jg.set_target_num_jobs_for_node( 5, 10 );

		TS_ASSERT_EQUALS( jg.get_num_nodes(), 5 );

		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 1 ),  1 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   1 ),  5 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 2 ),  6 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   2 ), 15 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 3 ), 26 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   3 ), 35 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 4 ), 16 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   4 ), 25 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 5 ), 36 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   5 ), 45 );

	}

	void test_genealogist_append_parents_and_n_replicate_jobs_for_node() {
		JobGenealogist jg;
		jg.set_num_nodes( 3 );
		jg.set_target_num_jobs_for_node( 1, 10 );
		jg.set_target_num_jobs_for_node( 2, 30 );
		jg.set_target_num_jobs_for_node( 3, 30 );
		jg.set_actual_njobs_for_node( 1, 10 );
		// soon! utility::vector1< core::Size > parents1 { 1, 5, 10 };
		utility::vector1< core::Size > parents1( 3 );
		parents1[ 1 ] = 1; parents1[ 2 ] = 5; parents1[ 3 ] = 10;
		jg.append_parents_and_n_replicate_jobs_for_node( 2, parents1, 10 );

		// soon! utility::vector1< core::Size > parents2 = { 12, 23 };
		utility::vector1< core::Size > parents2( 2 );
		parents2[ 1 ] = 12; parents2[ 2 ] = 23;
		jg.append_parents_and_n_replicate_jobs_for_node( 3, parents2, 10 );

		TS_ASSERT_EQUALS( jg.get_num_nodes(), 3 );

		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 1 ),  1 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   1 ), 10 );
		TS_ASSERT_EQUALS( jg.get_node_actual_range_begin( 1 ),  1 );
		TS_ASSERT_EQUALS( jg.get_node_actual_range_end(   1 ), 10 );

		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 2 ), 11 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   2 ), 40 );
		TS_ASSERT_EQUALS( jg.get_node_actual_range_begin( 2 ), 11 );
		TS_ASSERT_EQUALS( jg.get_node_actual_range_end(   2 ), 40 );

		TS_ASSERT_EQUALS( jg.get_node_target_range_begin( 3 ), 41 );
		TS_ASSERT_EQUALS( jg.get_node_target_range_end(   3 ), 70 );
		TS_ASSERT_EQUALS( jg.get_node_actual_range_begin( 3 ), 41 );
		TS_ASSERT_EQUALS( jg.get_node_actual_range_end(   3 ), 60 );

		for ( core::Size ii =  1; ii <= 10; ++ii ) { TS_ASSERT_EQUALS( jg.node_for_jobid( ii ), 1 );	}
		for ( core::Size ii = 11; ii <= 40; ++ii ) { TS_ASSERT_EQUALS( jg.node_for_jobid( ii ), 2 ); }
		for ( core::Size ii = 41; ii <= 60; ++ii ) { TS_ASSERT_EQUALS( jg.node_for_jobid( ii ), 3 ); }

		TS_ASSERT( ! jg.job_has_any_parents( 5 ) );
		TS_ASSERT(   jg.job_has_any_parents( 11 ) );
		TS_ASSERT(   jg.job_has_any_parents( 25 ) );
		TS_ASSERT(   jg.job_has_any_parents( 39 ) );
		TS_ASSERT(   jg.job_has_any_parents( 41 ) );
		TS_ASSERT(   jg.job_has_any_parents( 43 ) );
		TS_ASSERT(   jg.job_has_any_parents( 50 ) );
		TS_ASSERT(   jg.job_has_any_parents( 55 ) );

		TS_ASSERT_EQUALS( jg.get_parent_for_job( 11 ), 1 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 19 ), 1 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 20 ), 1 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 21 ), 5 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 24 ), 5 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 29 ), 5 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 30 ), 5 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 31 ), 10 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 32 ), 10 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 40 ), 10 );

		TS_ASSERT_EQUALS( jg.get_parent_for_job( 41 ), 12 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 49 ), 12 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 50 ), 12 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 51 ), 23 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 54 ), 23 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 59 ), 23 );
		TS_ASSERT_EQUALS( jg.get_parent_for_job( 60 ), 23 );

		TS_ASSERT_EQUALS( jg.get_parents_for_job( 19 ).size(),  1 );
		TS_ASSERT_EQUALS( jg.get_parents_for_job( 19 )[ 1 ],    1 );
		TS_ASSERT_EQUALS( jg.get_parents_for_job( 60 ).size(),  1 );
		TS_ASSERT_EQUALS( jg.get_parents_for_job( 60 )[ 1 ],   23 );

		// ok -- so now, let's start marking jobs as discarded
		std::list< core::Size > to_discard1 = { 2, 3, 4, 6, 7, 8, 9 };
		for ( auto ii : to_discard1 ) jg.note_job_discarded( ii );

		std::list< core::Size > res1 = jg.find_descendentless_jobs_backwards_from_node( 1 );
		TS_ASSERT( res1.empty() );

		std::list< core::Size > to_discard2 = {
			11, 13, 14, 15, 16, 17, 18, 19, 20,	21, 22, 24, 25, 26,
			27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 };
		for ( auto ii : to_discard2 ) jg.note_job_discarded( ii );

		std::list< core::Size > res2 = jg.find_descendentless_jobs_backwards_from_node( 2 );
		TS_ASSERT( ! res2.empty() );
		TS_ASSERT_EQUALS( res2.size(),   1 );
		TS_ASSERT_EQUALS( res2.front(), 10 );

		std::list< core::Size > to_discard3 = {
			10, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
		  51, 52, 53, 54, 56, 57, 58, 59, 60 };
		for ( auto ii : to_discard3 ) jg.note_job_discarded( ii );

		std::list< core::Size > res3a = jg.find_descendentless_jobs_backwards_from_node( 3 );
		TS_ASSERT( ! res3a.empty() );
		TS_ASSERT_EQUALS( res3a.size(),   1 );
		TS_ASSERT_EQUALS( res3a.front(), 12 );

		jg.note_job_discarded( 12 );

		// a second call, as we've discarded another batch of jobs, but start the search
		// from node 3
		std::list< core::Size > res3b = jg.find_descendentless_jobs_backwards_from_node( 3 );
		TS_ASSERT( ! res3b.empty() );
		TS_ASSERT_EQUALS( res3b.size(),   1 );
		TS_ASSERT_EQUALS( res3b.front(),  1 );
	}

	void test_get_next_job_for_node() {
		JobGenealogist jg;
		jg.set_num_nodes( 3 );
		jg.set_target_num_jobs_for_node( 1, 10 );
		jg.set_target_num_jobs_for_node( 2, 30 );
		jg.set_target_num_jobs_for_node( 3, 30 );

		jg.set_actual_njobs_for_node( 1, 10 );
		// soon! utility::vector1< core::Size > parents1 { 1, 5, 10 };
		utility::vector1< core::Size > parents1( 3 );
		parents1[ 1 ] = 1; parents1[ 2 ] = 5; parents1[ 3 ] = 10;
		jg.append_parents_and_n_replicate_jobs_for_node( 2, parents1, 10 );

		// soon! utility::vector1< core::Size > parents2 = { 12, 23 };
		utility::vector1< core::Size > parents2( 2 );
		parents2[ 1 ] = 12; parents2[ 2 ] = 23;
		jg.append_parents_and_n_replicate_jobs_for_node( 3, parents2, 10 );

		for ( core::Size ii = 1; ii <= 10; ++ii ) {
			core::Size next_job_id = jg.get_next_job_for_node( 1 );
			TS_ASSERT_EQUALS( next_job_id, ii );
		}
		TS_ASSERT_EQUALS( jg.get_next_job_for_node( 1 ), 0 );

		for ( core::Size ii = 11; ii <= 40; ++ii ) {
			core::Size next_job_id = jg.get_next_job_for_node( 2 );
			TS_ASSERT_EQUALS( next_job_id, ii );
		}
		TS_ASSERT_EQUALS( jg.get_next_job_for_node( 2 ), 0 );

		for ( core::Size ii = 41; ii <= 60; ++ii ) {
			core::Size next_job_id = jg.get_next_job_for_node( 3 );
			TS_ASSERT_EQUALS( next_job_id, ii );
		}
		TS_ASSERT_EQUALS( jg.get_next_job_for_node( 3 ), 0 );

	}

};
