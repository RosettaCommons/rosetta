// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/NodeManagers.cxxtest.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
//#include <test/core/init_util.hh>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/jd3/dag_node_managers/NodeManager.hh>
#include <protocols/jd3/dag_node_managers/SimpleNodeManager.hh>
#include <protocols/jd3/dag_node_managers/NodeManagerStorageMatrix.hh>

#include <basic/Tracer.hh>
static basic::Tracer TR("protocols.jd3.JD3NodeManagersTests");

using namespace protocols::jd3;
using namespace protocols::jd3::dag_node_managers;

class JD3NodeManagersTests : public CxxTest::TestSuite
{
public:

	JD3NodeManagersTests(){}

	void setUp() {
		//core_init();
		protocols_init();
	}

	void test_stop_early( ){

		SimpleNodeManager manager(
			0, // job offset
			4, // num jobs total
			3, // num results to keep
			7  // result threshold
		);

		//job 1 has 4 results
		//job 2 has 0 results
		//job 3 has 3 results
		//job 4 has 5 results

		manager.get_next_local_jobid();
		manager.note_job_completed( 1, 4 );
		manager.register_result( 1, 1, -50 );
		manager.register_result( 1, 2, -40 );
		manager.register_result( 1, 3, 3 );
		manager.register_result( 1, 4, 4 );

		manager.get_next_local_jobid();
		manager.note_job_completed( 2, 0 );

		manager.get_next_local_jobid();
		manager.note_job_completed( 3, 3 );
		manager.register_result( 3, 1, -30 );
		manager.register_result( 3, 2, 4 );
		manager.register_result( 3, 3, 4 );

		utility::vector1< ResultElements > const & elems = manager.results_to_keep();
		TS_ASSERT( manager.done_submitting() );

		manager.get_next_local_jobid();
		manager.note_job_completed( 4, 5 );
		manager.register_result( 4, 1, -60 );
		manager.register_result( 4, 2, -70 );
		manager.register_result( 4, 3, -80 );
		manager.register_result( 4, 4, -90 );
		manager.register_result( 4, 5, -99 );

		TS_ASSERT_EQUALS( elems.size(), 3 );

		int sum = 0;
		for ( ResultElements const & elem : elems ) {
			sum += core::Size( elem.score );
		}
		TS_ASSERT_EQUALS( sum, int( -50 - 40 - 30 ) );

		std::list< std::pair< core::Size, core::Size > > list;
		manager.append_job_results_that_should_be_discarded( list );
		TS_ASSERT_EQUALS( list.size(), 9 );
	}

	void single_node_manager_0_100_5( NodeManager & manager ){
		//NodeManager manager( 0, 100, 5 );

		TS_ASSERT_EQUALS( manager.job_offset(), 0 );
		TS_ASSERT_EQUALS( manager.num_jobs(), 100 );
		TS_ASSERT_EQUALS( manager.num_results_to_keep(), 5 );

		for ( core::Size ii = 1; ii <= 100; ++ii ) {
			TS_ASSERT_EQUALS( manager.num_jobs_submitted(), ii - 1 );
			TS_ASSERT( ! manager.done_submitting() );
			TS_ASSERT_EQUALS( ii, manager.get_next_local_jobid() );
			TS_ASSERT_EQUALS( manager.num_jobs_submitted(), ii );
			//TS_ASSERT( ! manager.all_results_are_in() );

			TS_ASSERT_EQUALS( manager.num_jobs_completed(), ii - 1 );
			manager.note_job_completed( ii, 1 );
			TS_ASSERT_EQUALS( manager.num_jobs_completed(), ii );
			manager.register_result( ii, 1, core::Real( 0 ) - ii, 1 );
		}

		TS_ASSERT( manager.all_results_are_in() );
		TS_ASSERT( manager.done_submitting() );

		TS_ASSERT( ! manager.all_waste_is_discarded() );
		std::list< std::pair< core::Size, core::Size > > dummy_list;
		manager.append_job_results_that_should_be_discarded( dummy_list );
		TS_ASSERT( manager.all_waste_is_discarded() );

		utility::vector1< ResultElements > const & results_to_keep = manager.results_to_keep();
		TS_ASSERT_EQUALS( results_to_keep.size(), 5 );
		TS_ASSERT_EQUALS( results_to_keep[ 1 ].score, -100 );
		TS_ASSERT_EQUALS( results_to_keep[ 2 ].score, -99 );
		TS_ASSERT_EQUALS( results_to_keep[ 3 ].score, -98 );
		TS_ASSERT_EQUALS( results_to_keep[ 4 ].score, -97 );
		TS_ASSERT_EQUALS( results_to_keep[ 5 ].score, -96 );
	}

	void test_single_node_manager_simple(){
		TR << "starting test_single_node_manager_simple()" << std::endl;
		SimpleNodeManager nm( 0, 100, 5 );
		single_node_manager_0_100_5( nm );

		utility::vector1< core::Size > dummy( 1, 100 );
		EvenlyPartitionedNodeManager snm ( 0, 100, 5, 1 );
		single_node_manager_0_100_5( snm );
	}

	void test_node_manager(){
		TR << "starting test_node_manager()" << std::endl;
		//core::Size const num_nodes = 3;

		utility::vector1< core::Size > num_jobs_for_node;
		num_jobs_for_node.push_back( 5 );
		num_jobs_for_node.push_back( 10 );
		num_jobs_for_node.push_back( 1 );

		utility::vector1< core::Size > num_results_to_keep_for_node;
		num_results_to_keep_for_node.push_back( 10 );
		num_results_to_keep_for_node.push_back( 5 );
		num_results_to_keep_for_node.push_back( 1 );

		utility::vector1< NodeManagerOP > managers;
		managers.reserve( 3 );
		managers.push_back( NodeManagerOP( new SimpleNodeManager( 0, num_jobs_for_node[ 1 ], num_results_to_keep_for_node[ 1 ] ) ) );
		managers.push_back( NodeManagerOP( new SimpleNodeManager( num_jobs_for_node[ 1 ], num_jobs_for_node[ 2 ], num_results_to_keep_for_node[ 2 ] ) ) );
		managers.push_back( NodeManagerOP( new SimpleNodeManager( num_jobs_for_node[ 1 ] + num_jobs_for_node[ 2 ], num_jobs_for_node[ 3 ], num_results_to_keep_for_node[ 3 ] ) ) );

		//Node 1 focuses on having multiple results per job
		utility::vector1< utility::vector1< core::Real > > score_for_result_for_local_job_id_node1;
		score_for_result_for_local_job_id_node1.resize( 5 );

		score_for_result_for_local_job_id_node1[ 1 ].push_back( -0.2 );
		score_for_result_for_local_job_id_node1[ 1 ].push_back( -0.4 );//KEEP 10
		score_for_result_for_local_job_id_node1[ 1 ].push_back( -0.1 );

		score_for_result_for_local_job_id_node1[ 2 ].push_back( -2.1 );//KEEP 3

		score_for_result_for_local_job_id_node1[ 3 ].push_back( -4 );//KEEP 1
		score_for_result_for_local_job_id_node1[ 3 ].push_back( 1.2 );
		score_for_result_for_local_job_id_node1[ 3 ].push_back( 0 );
		score_for_result_for_local_job_id_node1[ 3 ].push_back( -0.6 );//KEEP 8
		score_for_result_for_local_job_id_node1[ 3 ].push_back( -1.0 );//KEEP 4

		score_for_result_for_local_job_id_node1[ 4 ].push_back( -0.3 );
		score_for_result_for_local_job_id_node1[ 4 ].push_back( -0.5 );//KEEP 9
		score_for_result_for_local_job_id_node1[ 4 ].push_back( -0.7 );//KEEP 7

		score_for_result_for_local_job_id_node1[ 5 ].push_back( -0.9 );//KEEP 5/6
		score_for_result_for_local_job_id_node1[ 5 ].push_back( -0.9 );//KEEP 5/6
		score_for_result_for_local_job_id_node1[ 5 ].push_back( -3.3 );//KEEP 2

		{
			NodeManager & manager = * managers[ 1 ];
			for ( core::Size local_job_id = 1; local_job_id <= 5; ++local_job_id ) {
				TS_ASSERT_EQUALS( manager.num_jobs_submitted(), local_job_id - 1 );
				TS_ASSERT( ! manager.done_submitting() );
				TS_ASSERT_EQUALS( local_job_id, manager.get_next_local_jobid() );
				TS_ASSERT_EQUALS( manager.num_jobs_submitted(), local_job_id );
				//TS_ASSERT( ! manager.all_results_are_in() );

				core::Size const num_results = score_for_result_for_local_job_id_node1[ local_job_id ].size();
				manager.note_job_completed( local_job_id, num_results );
				for ( core::Size result_id = 1; result_id <= num_results; ++result_id ) {
					manager.register_result( local_job_id, result_id, score_for_result_for_local_job_id_node1[ local_job_id ][ result_id ], 1 );
				}

			}

			TS_ASSERT( manager.done_submitting() );
			TS_ASSERT( manager.all_results_are_in() );

			TS_ASSERT( ! manager.all_waste_is_discarded() );
			std::list< std::pair< core::Size, core::Size > > dummy_list;
			manager.append_job_results_that_should_be_discarded( dummy_list );
			TS_ASSERT_EQUALS( dummy_list.size(), 5 );
			TS_ASSERT( manager.all_waste_is_discarded() );

			utility::vector1< ResultElements > const & results_to_keep = manager.results_to_keep();
			TS_ASSERT_EQUALS( results_to_keep.size(), 10 );
			TS_ASSERT_EQUALS( results_to_keep[ 1 ].score, -4 );
			TS_ASSERT_EQUALS( results_to_keep[ 2 ].score, -3.3 );
			TS_ASSERT_EQUALS( results_to_keep[ 3 ].score, -2.1 );
			TS_ASSERT_EQUALS( results_to_keep[ 4 ].score, -1.0 );
			TS_ASSERT_EQUALS( results_to_keep[ 5 ].score, -0.9 );
			TS_ASSERT_EQUALS( results_to_keep[ 6 ].score, -0.9 );
			TS_ASSERT_EQUALS( results_to_keep[ 7 ].score, -0.7 );
			TS_ASSERT_EQUALS( results_to_keep[ 8 ].score, -0.6 );
			TS_ASSERT_EQUALS( results_to_keep[ 9 ].score, -0.5 );
			TS_ASSERT_EQUALS( results_to_keep[ 10 ].score, -0.4 );

		}

		//Node 2 focuses on having some runs with 0 results
		utility::vector1< utility::vector1< core::Real > > score_for_result_for_local_job_id_node2;
		score_for_result_for_local_job_id_node2.resize( 10 );

		score_for_result_for_local_job_id_node2[ 1 ].push_back( -1 );//5

		score_for_result_for_local_job_id_node2[ 3 ].push_back( -3 );//3
		score_for_result_for_local_job_id_node2[ 4 ].push_back( -4 );//2
		score_for_result_for_local_job_id_node2[ 5 ].push_back( -2 );//4
		score_for_result_for_local_job_id_node2[ 6 ].push_back( -0.1 );

		score_for_result_for_local_job_id_node2[ 8 ].push_back( -0.31 );
		score_for_result_for_local_job_id_node2[ 9 ].push_back( -10 );//1

		{
			NodeManager & manager = * managers[ 2 ];
			for ( core::Size local_job_id = 1; local_job_id <= 10; ++local_job_id ) {
				core::Size global_job_id = num_jobs_for_node[ 1 ] + local_job_id;

				TS_ASSERT_EQUALS( manager.num_jobs_submitted(), local_job_id - 1 );
				TS_ASSERT( ! manager.done_submitting() );
				TS_ASSERT_EQUALS( local_job_id, manager.get_next_local_jobid() );
				TS_ASSERT_EQUALS( manager.num_jobs_submitted(), local_job_id );
				//TS_ASSERT( ! manager.all_results_are_in() );

				core::Size const num_results = score_for_result_for_local_job_id_node2[ local_job_id ].size();

				manager.note_job_completed( global_job_id, num_results );
				for ( core::Size result_id = 1; result_id <= num_results; ++result_id ) {
					manager.register_result(
						global_job_id,
						result_id,
						score_for_result_for_local_job_id_node2[ local_job_id ][ result_id ],
						1
					);
				}

			}

			TS_ASSERT( manager.done_submitting() );
			TS_ASSERT( manager.all_results_are_in() );

			TS_ASSERT( ! manager.all_waste_is_discarded() );
			std::list< std::pair< core::Size, core::Size > > dummy_list;
			manager.append_job_results_that_should_be_discarded( dummy_list );
			TS_ASSERT_EQUALS( dummy_list.size(), 2 );
			TS_ASSERT( manager.all_waste_is_discarded() );

			utility::vector1< ResultElements > const & results_to_keep = manager.results_to_keep();
			TS_ASSERT_EQUALS( results_to_keep.size(), 5 );
			TS_ASSERT_EQUALS( results_to_keep[ 1 ].score, -10 );
			TS_ASSERT_EQUALS( results_to_keep[ 2 ].score, -4 );
			TS_ASSERT_EQUALS( results_to_keep[ 3 ].score, -3 );
			TS_ASSERT_EQUALS( results_to_keep[ 4 ].score, -2 );
			TS_ASSERT_EQUALS( results_to_keep[ 5 ].score, -1 );

		}


		//Node 3 has rare case of 1 job and 1 result to keep
		utility::vector1< utility::vector1< core::Real > > score_for_result_for_local_job_id_node3;
		score_for_result_for_local_job_id_node3.resize( 1 );

		score_for_result_for_local_job_id_node3[ 1 ].push_back( -4 );
		score_for_result_for_local_job_id_node3[ 1 ].push_back( -3 );
		score_for_result_for_local_job_id_node3[ 1 ].push_back( -2 );
		score_for_result_for_local_job_id_node3[ 1 ].push_back( -1 );
		score_for_result_for_local_job_id_node3[ 1 ].push_back( -1 );

		{
			NodeManager & manager = * managers[ 3 ];
			for ( core::Size local_job_id = 1; local_job_id <= 1; ++local_job_id ) {
				core::Size global_job_id = num_jobs_for_node[ 1 ] + num_jobs_for_node[ 2 ] + local_job_id;

				TS_ASSERT_EQUALS( manager.num_jobs_submitted(), local_job_id - 1 );
				TS_ASSERT( ! manager.done_submitting() );
				TS_ASSERT_EQUALS( local_job_id, manager.get_next_local_jobid() );
				TS_ASSERT_EQUALS( manager.num_jobs_submitted(), local_job_id );
				//TS_ASSERT( ! manager.all_results_are_in() );

				core::Size const num_results = score_for_result_for_local_job_id_node3[ local_job_id ].size();

				manager.note_job_completed( global_job_id, num_results );
				for ( core::Size result_id = 1; result_id <= num_results; ++result_id ) {
					manager.register_result( global_job_id, result_id, score_for_result_for_local_job_id_node3[ local_job_id ][ result_id ], 1 );
				}

			}

			TS_ASSERT( manager.done_submitting() );
			TS_ASSERT( manager.all_results_are_in() );

			TS_ASSERT( ! manager.all_waste_is_discarded() );
			std::list< std::pair< core::Size, core::Size > > dummy_list;
			manager.append_job_results_that_should_be_discarded( dummy_list );
			TS_ASSERT_EQUALS( dummy_list.size(), 4 );
			TS_ASSERT( manager.all_waste_is_discarded() );

			utility::vector1< ResultElements > const & results_to_keep = manager.results_to_keep();
			TS_ASSERT_EQUALS( results_to_keep.size(), 1 );
			TS_ASSERT_EQUALS( results_to_keep[ 1 ].score, -4 );

		}
	}

	struct result{
		result( core::Size gji, core::Size lri, core::Real sc, core::Size ipi ){
			global_job_id = gji;
			local_result_id = lri;
			score = sc;
			input_pose_id = ipi;
		}

		core::Size global_job_id;
		core::Size local_result_id;
		core::Real score;
		core::Size input_pose_id;
	};

	void test_segregated_node_manager(){
		TR << "starting test_segregated_node_manager()" << std::endl;

		core::Size const job_offset = 100;
		core::Size const num_jobs_total = 6;
		core::Size const num_results_to_keep = 5;
		core::Size const num_segregations = 2;

		EvenlyPartitionedNodeManager manager( job_offset, num_jobs_total, num_results_to_keep, num_segregations );
		for ( core::Size global_job_id = 101; global_job_id <= 110; ++global_job_id ) {
			manager.note_job_completed( global_job_id, 2 );
		}

		std::list< result > results;
		results.push_back( result( 101, 1, -10, 1 ) );
		results.push_back( result( 101, 2, 112, 1 ) );
		results.push_back( result( 102, 1, -50, 1 ) );//winner
		results.push_back( result( 102, 2, -30, 1 ) );
		results.push_back( result( 103, 1, -20, 1 ) );
		results.push_back( result( 103, 2, -11, 1 ) );
		results.push_back( result( 104, 1, -80, 1 ) );//winner
		results.push_back( result( 104, 2, -35, 1 ) );
		results.push_back( result( 105, 1, -25, 1 ) );
		results.push_back( result( 105, 2, -54, 1 ) );//winner

		results.push_back( result( 106, 1, 90, 2 ) );
		results.push_back( result( 106, 2, 112, 2 ) );
		results.push_back( result( 107, 1, 50, 2 ) );
		results.push_back( result( 107, 2, 30, 2 ) );
		results.push_back( result( 108, 1, 20, 2 ) );//winner
		results.push_back( result( 108, 2, 10, 2 ) );//winner
		results.push_back( result( 109, 1, 80, 2 ) );
		results.push_back( result( 109, 2, 35, 2 ) );
		results.push_back( result( 110, 1, 25, 2 ) );
		results.push_back( result( 110, 2, 54, 2 ) );

		for ( result ii : results ) {
			manager.register_result( ii.global_job_id, ii.local_result_id, ii.score, ii.input_pose_id );
		}

		utility::vector1< ResultElements > results_to_keep = manager.results_to_keep();
		TS_ASSERT_EQUALS( results_to_keep.size(), num_results_to_keep );
		int sum = 0;
		for ( ResultElements & elem : results_to_keep ) {
			sum += int( elem.score );
		}
		TS_ASSERT_EQUALS( sum, 20 + 10 - 50 - 80 - 54);
	}

	void test_simple1(){
		unsigned int const job_offset = 0;
		//num_jobs_total = 100
		//num_results_to_keep = 3
		unsigned int const result_threshold = 10;

		SimpleNodeManager snm( 0, 100, 3, result_threshold );
		TS_ASSERT_EQUALS( snm.job_offset(), job_offset );

		for ( unsigned int i = 1; i <= 90; ++i ) {
			TS_ASSERT_EQUALS( snm.get_next_local_jobid(), i );
			TS_ASSERT_EQUALS( snm.num_jobs_submitted(), i );
		}

		for ( unsigned int i = 1; i <= 90; ++i ) {
			snm.note_job_completed( i, 1 );
			unsigned int const offset_i = (i + 10) % 90;
			snm.register_result( i, 1, -1.0 * offset_i );
			TS_ASSERT_EQUALS( snm.done_submitting(), i >= result_threshold );
			TS_ASSERT_EQUALS( snm.jobs_are_still_running(), i != 90 );
		}

		//    i    score
		//    1    -11
		//    2    -12
		//   25    -35
		//   79    -89
		//   80    0
		//   81    -1
		//
		//Best 3 scores are -89, -88, -87

		std::list< std::pair< core::Size, core::Size > > list;
		snm.append_job_results_that_should_be_discarded( list );
		TS_ASSERT_EQUALS( list.size(), 90 - 3 );

		utility::vector1< ResultElements > const & results = snm.results_to_keep();
		TS_ASSERT_DELTA( results[ 1 ].score, -89, 0.1 );
		TS_ASSERT_DELTA( results[ 2 ].score, -88, 0.1 );
		TS_ASSERT_DELTA( results[ 3 ].score, -87, 0.1 );
	}

	void test_storage_matric_depth_first(){
		/*
		Structure:

		partition
		1          -2  -1   x
		2          -4  -3  -1

		*/

		using namespace protocols::jd3::dag_node_managers;

		utility::vector1< core::Size > n_results_to_keep_for_partition;
		n_results_to_keep_for_partition.push_back( 3 );
		n_results_to_keep_for_partition.push_back( 3 );
		bool depth_first = true;

		NodeManagerStorageMatrix matrix( n_results_to_keep_for_partition, depth_first );

		matrix.insert( 1, ResultElements( 1, 1, -1 ) );
		matrix.insert( 1, ResultElements( 1, 1, -2 ) );

		matrix.insert( 2, ResultElements( 1, 1, 3 ) );
		matrix.insert( 2, ResultElements( 1, 1, -1 ) );
		matrix.insert( 2, ResultElements( 1, 1, -4 ) );
		matrix.insert( 2, ResultElements( 1, 1, 7 ) );
		matrix.insert( 2, ResultElements( 1, 1, -3 ) );
		matrix.insert( 2, ResultElements( 1, 1, 5 ) );

		TS_ASSERT_EQUALS( matrix.get_nth_element( 1 ).score, -2 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 2 ).score, -1 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 3 ).score, 0 );

		TS_ASSERT_EQUALS( matrix.get_nth_element( 4 ).score, -4 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 5 ).score, -3 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 6 ).score, -1 );
	}

	void test_storage_matric_breadth_first(){
		/*
		Structure:

		partition
		1          -2  -1   x
		2          -4  -3  -1

		*/

		using namespace protocols::jd3::dag_node_managers;

		utility::vector1< core::Size > n_results_to_keep_for_partition;
		n_results_to_keep_for_partition.push_back( 4 );
		n_results_to_keep_for_partition.push_back( 4 );
		bool depth_first = false;

		NodeManagerStorageMatrix matrix( n_results_to_keep_for_partition, depth_first );

		matrix.insert( 1, ResultElements( 1, 1, -1 ) );
		matrix.insert( 1, ResultElements( 1, 1, -2 ) );

		matrix.insert( 2, ResultElements( 1, 1, -1 ) );
		matrix.insert( 2, ResultElements( 1, 1, -4 ) );
		matrix.insert( 2, ResultElements( 1, 1, -3 ) );

		TS_ASSERT_EQUALS( matrix.get_nth_element( 1 ).score, -2 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 3 ).score, -1 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 5 ).score, 0 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 2 ).score, -4 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 4 ).score, -3 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 6 ).score, -1 );

		matrix.insert( 1, ResultElements( 1, 1, -6 ) );

		matrix.insert( 2, ResultElements( 1, 1, -6 ) );
		matrix.insert( 2, ResultElements( 1, 1, -7 ) );

		TS_ASSERT_EQUALS( matrix.get_nth_element( 7 ).score, -6 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 8 ).score, -7 );
	}

	void test_storage_matric_breadth_first2(){
		/*
		Structure:

		partition
		1          -2  -1   x
		2          -4  -3  -1

		*/

		using namespace protocols::jd3::dag_node_managers;

		utility::vector1< core::Size > n_results_to_keep_for_partition;
		n_results_to_keep_for_partition.push_back( 3 );
		n_results_to_keep_for_partition.push_back( 4 );
		bool depth_first = false;

		NodeManagerStorageMatrix matrix( n_results_to_keep_for_partition, depth_first );

		matrix.insert( 1, ResultElements( 1, 1, -1 ) );
		matrix.insert( 1, ResultElements( 1, 1, -2 ) );

		matrix.insert( 2, ResultElements( 1, 1, -1 ) );
		matrix.insert( 2, ResultElements( 1, 1, -4 ) );
		matrix.insert( 2, ResultElements( 1, 1, -3 ) );

		TS_ASSERT_EQUALS( matrix.get_nth_element( 1 ).score, -2 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 3 ).score, -1 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 5 ).score, 0 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 2 ).score, -4 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 4 ).score, -3 );
		TS_ASSERT_EQUALS( matrix.get_nth_element( 6 ).score, -1 );

		matrix.insert( 1, ResultElements( 1, 1, -6 ) );

		matrix.insert( 2, ResultElements( 1, 1, -6 ) );
		matrix.insert( 2, ResultElements( 1, 1, -7 ) );

		TS_ASSERT_EQUALS( matrix.get_nth_element( 7 ).score, -7 );
	}

	void test_tokens(){
		using namespace protocols::jd3::dag_node_managers;
		utility::vector1< core::Real > scores_for_token1;
		utility::vector1< core::Real > scores_for_token2;
		utility::vector1< core::Real > scores_for_token3;

		scores_for_token1.push_back( -1.0 );
		scores_for_token1.push_back( -2.0 );
		scores_for_token1.push_back( -3.0 );
		scores_for_token1.push_back( -4.0 );

		scores_for_token2.push_back( 1.0 );
		scores_for_token2.push_back( 2.0 );
		scores_for_token2.push_back( 3.0 );
		scores_for_token2.push_back( 4.0 );

		scores_for_token3.push_back( 4.1 );
		scores_for_token3.push_back( 4.4 );
		scores_for_token3.push_back( 2.2 );
		scores_for_token3.push_back( -3.5 );

		utility::vector1< core::Size > const n_results_to_keep_for_partition( 1, 5 );
		bool depth_first = false;
		core::Size const max_num_results_with_same_token_per_partition = 2;

		NodeManagerStorageMatrix matrix( n_results_to_keep_for_partition, depth_first );
		matrix.set_max_num_results_with_same_token_per_partition( max_num_results_with_same_token_per_partition );

		TS_ASSERT_EQUALS( matrix.num_partitions(), 1 );

		//Results should look like:
		//Score: -4.0 -3.5 -3.0 1.0 2.0
		//Rank:  1    2    3    4   5
		//Token: 1    3    1    2   2

		for ( int i=1; i<5; ++i ) {
			matrix.insert( 1, ResultElements( 1, i, scores_for_token1[ i ], 1 ) );
			matrix.insert( 1, ResultElements( 2, i, scores_for_token2[ i ], 2 ) );
			matrix.insert( 1, ResultElements( 3, i, scores_for_token3[ i ], 3 ) );
		}

		utility::vector1< ResultElements > linear_vector_of_results = matrix.linear_vector_of_results();

		TS_ASSERT_EQUALS( linear_vector_of_results[ 1 ].token, 1 );
		TS_ASSERT_EQUALS( linear_vector_of_results[ 2 ].token, 3 );
		TS_ASSERT_EQUALS( linear_vector_of_results[ 3 ].token, 1 );
		TS_ASSERT_EQUALS( linear_vector_of_results[ 4 ].token, 2 );
		TS_ASSERT_EQUALS( linear_vector_of_results[ 5 ].token, 2 );
	}

private:

};
