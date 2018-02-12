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
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/jd3/JobGenealogist.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>
#include <basic/Tracer.hh>

static basic::Tracer TR("test.protocols.jd3.JobGenealogistTests");

using namespace protocols::jd3;

///@brief dummy class that gives us access to protected methods
class JobGenealogist_ : public JobGenealogist {

public:

	JobGenealogist_(
		core::Size num_job_dag_nodes,
		core::Size num_input_sources
	) :
		JobGenealogist( num_job_dag_nodes, num_input_sources ){}

	JGJobNode const * _get_job_node( core::Size job_dag_node, core::Size global_job_id ) const {
		return JobGenealogist::get_job_node( job_dag_node, global_job_id );
	}

	JGResultNode const * _get_result_node( core::Size node, core::Size global_job_id, core::Size result_id ) const{
		return JobGenealogist::get_result_node( node, global_job_id, result_id );
	}

	JGJobNode * _get_job_node( core::Size job_dag_node, core::Size global_job_id ){
		return JobGenealogist::get_job_node( job_dag_node, global_job_id );
	}

	JGResultNode * _get_result_node( core::Size node, core::Size global_job_id, core::Size result_id ) {
		return JobGenealogist::get_result_node( node, global_job_id, result_id );
	}

};

class JobGenealogistTests : public CxxTest::TestSuite
{
public:

	void setUp() {
		protocols_init();
	}

	void test_garbage_collection1(){
		TR << "beginning test_garbage_collection1" << std::endl;
		inner_test_garbage_collection( true );
	}

	void test_garbage_collection2(){
		TR << "beginning test_garbage_collection2" << std::endl;
		inner_test_garbage_collection( false );
	}

	void inner_test_garbage_collection(
		bool optionA
	){

		/* TREE:

		J1    J2    J3    J4
		|     |     |
		R1    R2    R3
		|     |
		J5    J6
		|
		R4
		*/


		JobGenealogist_ jobgen( 2, 1 );
		auto job1 = jobgen.register_new_job( 1, 1, 1 );
		TS_ASSERT( job1 );
		auto job2 = jobgen.register_new_job( 1, 2, 1 );
		TS_ASSERT( job2 );
		auto job3 = jobgen.register_new_job( 1, 3, 1 );
		TS_ASSERT( job3 );
		auto job4 = jobgen.register_new_job( 1, 4, 1 );
		TS_ASSERT( job4 );

		jobgen.note_job_completed( job1, 1 );
		jobgen.note_job_completed( job2, 1 );
		jobgen.note_job_completed( job3, 1 );
		jobgen.note_job_completed( job4, 0 );

		TS_ASSERT_EQUALS( job1->children().size(), 1 );
		TS_ASSERT_EQUALS( job2->children().size(), 1 );
		TS_ASSERT_EQUALS( job3->children().size(), 1 );
		TS_ASSERT_EQUALS( job4->children().size(), 0 );

		auto result1 = job1->children()[ 1 ];
		auto result2 = job2->children()[ 1 ];
		auto result3 = job3->children()[ 1 ];

		TS_ASSERT( result1 );
		TS_ASSERT( result2 );
		TS_ASSERT( result3 );

		auto job5 = jobgen.register_new_job( 2, 5, result1 );
		auto job6 = jobgen.register_new_job( 2, 6, result2 );

		TS_ASSERT( job5 );
		TS_ASSERT( job6 );

		jobgen.note_job_completed( job5, 1 );
		jobgen.note_job_completed( job6, 0 );

		std::list< std::pair< core::Size, core::Size > > container_for_discarded_result_ids;
		jobgen.garbage_collection( 1, optionA, container_for_discarded_result_ids );

		TS_ASSERT_EQUALS( job1->children().size(), 1 );
		TS_ASSERT_EQUALS( job2->children().size(), (optionA ? 0 : 1) );
		TS_ASSERT_EQUALS( job3->children().size(), 0 );
		TS_ASSERT_EQUALS( job4->children().size(), 0 );
		if ( optionA ) {
			TS_ASSERT_EQUALS( container_for_discarded_result_ids.size(), 2 );

			auto iter = std::find(
				container_for_discarded_result_ids.begin(),
				container_for_discarded_result_ids.end(),
				protocols::jd3::JobResultID{ 2, 1 }
			);
			TS_ASSERT( iter != container_for_discarded_result_ids.end() );

			iter = std::find(
				container_for_discarded_result_ids.begin(),
				container_for_discarded_result_ids.end(),
				protocols::jd3::JobResultID{ 3, 1 }
			);
			TS_ASSERT( iter != container_for_discarded_result_ids.end() );

		} else {
			TS_ASSERT_EQUALS( container_for_discarded_result_ids.size(), 1 );

			auto iter = std::find(
				container_for_discarded_result_ids.begin(),
				container_for_discarded_result_ids.end(),
				protocols::jd3::JobResultID{ 3, 1 }
			);
			TS_ASSERT( iter != container_for_discarded_result_ids.end() );
		}

	}

	void test_garbage_collection_with_multiple_results1(){
		TR << "beginning test_garbage_collection_with_multiple_results1" << std::endl;
		inner_test_garbage_collection_with_multiple_results( true );
	}

	void test_garbage_collection_with_multiple_results2(){
		TR << "beginning test_garbage_collection_with_multiple_results2" << std::endl;
		inner_test_garbage_collection_with_multiple_results( false );
	}

	void inner_test_garbage_collection_with_multiple_results(
		bool delete_downstream_job_if_it_has_no_results
	){
		// Tree. Results with captial R should be deleted and lowercase r should be kept.
		// ? is for results that should be deleted if and only if delete_downstream_job_if_it_has_no_results == true

		// J1     J2     J3     J4      NODE 1
		// |      |      |      |
		// | \    | \    | \    | \_
		// |  |   |  |   |  |   |  |
		// r1 ?2  r1 R2  R1 ?2  R1 R2
		// |  |   |         |
		// J5 J6  J7        J8          NODE 2
		// |      | \_
		// r1     r1 r2

		JobGenealogist_ jobgen( 2, 1 );
		auto job1 = jobgen.register_new_job( 1, 1, 1 );
		TS_ASSERT( job1 );
		auto job2 = jobgen.register_new_job( 1, 2, 1 );
		TS_ASSERT( job2 );
		auto job3 = jobgen.register_new_job( 1, 3, 1 );
		TS_ASSERT( job3 );
		auto job4 = jobgen.register_new_job( 1, 4, 1 );
		TS_ASSERT( job4 );

		jobgen.note_job_completed( job1, 2 );
		jobgen.note_job_completed( job2, 2 );
		jobgen.note_job_completed( job3, 2 );
		jobgen.note_job_completed( job4, 2 );

		TS_ASSERT_EQUALS( job1->children().size(), 2 );
		TS_ASSERT_EQUALS( job2->children().size(), 2 );
		TS_ASSERT_EQUALS( job3->children().size(), 2 );
		TS_ASSERT_EQUALS( job4->children().size(), 2 );

		auto job5 = jobgen.register_new_job( 2, 5, 1, 1, 1 );//J5
		auto job6 = jobgen.register_new_job( 2, 6, 1, 1, 2 );//J6
		auto job7 = jobgen.register_new_job( 2, 7, 1, 2, 1 );//J7
		auto job8 = jobgen.register_new_job( 2, 8, 1, 3, 2 );//J8

		TS_ASSERT( job5 );
		TS_ASSERT( job6 );
		TS_ASSERT( job7 );
		TS_ASSERT( job8 );

		TS_ASSERT_EQUALS( job1->global_job_id(), 1 );
		TS_ASSERT_EQUALS( job2->global_job_id(), 2 );
		TS_ASSERT_EQUALS( job3->global_job_id(), 3 );
		TS_ASSERT_EQUALS( job4->global_job_id(), 4 );
		TS_ASSERT_EQUALS( job5->global_job_id(), 5 );
		TS_ASSERT_EQUALS( job6->global_job_id(), 6 );
		TS_ASSERT_EQUALS( job7->global_job_id(), 7 );
		TS_ASSERT_EQUALS( job8->global_job_id(), 8 );

		jobgen.note_job_completed( job5, 1 );
		jobgen.note_job_completed( job6, 0 );
		jobgen.note_job_completed( job7, 2 );
		jobgen.note_job_completed( job8, 0 );

		TS_ASSERT_EQUALS( job5->children().size(), 1 );
		TS_ASSERT_EQUALS( job6->children().size(), 0 );
		TS_ASSERT_EQUALS( job7->children().size(), 2 );
		TS_ASSERT_EQUALS( job8->children().size(), 0 );

		//It may seem silly that we are checking this again, but there used to be a bug that occurred since the last time we checked these
		//Figured I would leave them here. Who knows? It could happen again
		TS_ASSERT_EQUALS( job1->global_job_id(), 1 );
		TS_ASSERT_EQUALS( job2->global_job_id(), 2 );
		TS_ASSERT_EQUALS( job3->global_job_id(), 3 );
		TS_ASSERT_EQUALS( job4->global_job_id(), 4 );
		TS_ASSERT_EQUALS( job5->global_job_id(), 5 );
		TS_ASSERT_EQUALS( job6->global_job_id(), 6 );
		TS_ASSERT_EQUALS( job7->global_job_id(), 7 );
		TS_ASSERT_EQUALS( job8->global_job_id(), 8 );

		std::list< std::pair< core::Size, core::Size > > container_for_discarded_result_ids;
		jobgen.garbage_collection( 1, delete_downstream_job_if_it_has_no_results, container_for_discarded_result_ids );

		TS_ASSERT_EQUALS( job1->global_job_id(), 1 );
		TS_ASSERT_EQUALS( job2->global_job_id(), 2 );
		TS_ASSERT_EQUALS( job5->global_job_id(), 5 );
		TS_ASSERT_EQUALS( job7->global_job_id(), 7 );

		if ( ! delete_downstream_job_if_it_has_no_results ) {
			TS_ASSERT_EQUALS( job3->global_job_id(), 3 );
			TS_ASSERT_EQUALS( job6->global_job_id(), 6 );
			TS_ASSERT_EQUALS( job8->global_job_id(), 8 );
		}

		TS_ASSERT_EQUALS( container_for_discarded_result_ids.size(),
			( delete_downstream_job_if_it_has_no_results ? 6 : 4 ) );

		auto iter = std::find(
			container_for_discarded_result_ids.begin(),
			container_for_discarded_result_ids.end(),
			protocols::jd3::JobResultID{ 2, 2 }
		);
		TS_ASSERT( iter != container_for_discarded_result_ids.end() );

		iter = std::find(
			container_for_discarded_result_ids.begin(),
			container_for_discarded_result_ids.end(),
			protocols::jd3::JobResultID{ 3, 1 }
		);
		TS_ASSERT( iter != container_for_discarded_result_ids.end() );

		iter = std::find(
			container_for_discarded_result_ids.begin(),
			container_for_discarded_result_ids.end(),
			protocols::jd3::JobResultID{ 4, 1 }
		);
		TS_ASSERT( iter != container_for_discarded_result_ids.end() );

		iter = std::find(
			container_for_discarded_result_ids.begin(),
			container_for_discarded_result_ids.end(),
			protocols::jd3::JobResultID{ 4, 2 }
		);
		TS_ASSERT( iter != container_for_discarded_result_ids.end() );

		if ( delete_downstream_job_if_it_has_no_results ) {
			iter = std::find(
				container_for_discarded_result_ids.begin(),
				container_for_discarded_result_ids.end(),
				protocols::jd3::JobResultID{ 1, 2 }
			);
			TS_ASSERT( iter != container_for_discarded_result_ids.end() );

			iter = std::find(
				container_for_discarded_result_ids.begin(),
				container_for_discarded_result_ids.end(),
				protocols::jd3::JobResultID{ 3, 2 }
			);
			TS_ASSERT( iter != container_for_discarded_result_ids.end() );
		}

	}


	void test_multiple_parents(){
		TR << "beginning test_multiple_parents( true )" << std::endl;
		inner_test_multiple_parents( true );

		TR << "beginning test_multiple_parents( false )" << std::endl;
		inner_test_multiple_parents( false );
	}

	void inner_test_multiple_parents( bool optionA ) {
		JobGenealogist_ jobgen( 2, 2 );
		JobGenealogist_ const * const_jobgen = & jobgen;

		auto job1 = jobgen.register_new_job( 1, 1, 1 );
		TS_ASSERT( job1 );
		auto job2 = jobgen.register_new_job( 1, 2, 2 );
		TS_ASSERT( job2 );

		jobgen.note_job_completed( 1, 1, 1 );
		jobgen.note_job_completed( 1, 2, 1 );

		auto result1 = jobgen._get_result_node( 1, 1, 1 );
		auto result2 = jobgen._get_result_node( 1, 2, 1 );

		TS_ASSERT_EQUALS( const_jobgen->_get_result_node(1,1,1), result1 );

		decltype( job1 ) job3;
		if ( optionA ) {
			job3 = jobgen.register_new_job(
				2, //node_id
				3, //global_job_id
				1, //node_id_of_parent
				1, //global_job_id_of_parent
				1  //result_id_of_parent
			);
			jobgen.note_job_completed( 2, 3, 1 );
			TS_ASSERT_EQUALS( job3->parents().size(), 1 );
			TS_ASSERT_EQUALS( job3->parents()[ 1 ], result1 );
			TS_ASSERT_EQUALS( job3->input_source_id(), 1 );

			job3->add_parent( result2, true );

		} else {
			utility::vector1< JGResultNode * > parents;
			parents.push_back( result1 );
			parents.push_back( result2 );
			job3 = jobgen.register_new_job( 2, 3, parents );
			jobgen.note_job_completed( 2, 3, 1 );
		}

		TS_ASSERT_EQUALS( job3->parents().size(), 2 );
		TS_ASSERT_EQUALS( job3->parents()[ 1 ], result1 );
		TS_ASSERT_EQUALS( job3->parents()[ 2 ], result2 );
		TS_ASSERT_EQUALS( job3->input_source_id(), 1 );

		job3->remove_parent( result1, true );
		TS_ASSERT_EQUALS( job3->parents().size(), 1 );
		TS_ASSERT_EQUALS( job3->parents()[ 1 ], result2 );
		TS_ASSERT_EQUALS( job3->input_source_id(), 2 );

		//Check to make sure job1 has no link to job3
		std::list< std::pair< core::Size, core::Size > > container_for_discarded_result_ids;
		jobgen.garbage_collection( 1, true, container_for_discarded_result_ids );
		TS_ASSERT_EQUALS( container_for_discarded_result_ids.size(), 1 );
		std::pair< core::Size, core::Size > const & trash = * container_for_discarded_result_ids.begin();
		TS_ASSERT_EQUALS( trash.first, 1 );
		TS_ASSERT_EQUALS( trash.second, 1 );
	}

	void test_job_genealogist2_no_results_for_node1(){
		TR << "beginning test_job_genealogist2_no_results_for_node1()" << std::endl;

		unsigned int const num_job_dag_nodes = 3;
		unsigned int const num_input_sources = 2;
		JobGenealogist_ jobgen( num_job_dag_nodes, num_input_sources );

		unsigned int const num_jobs_node1 = 10;
		unsigned int const half_num_jobs_node1 = num_jobs_node1/2;
		for ( unsigned int i = 1; i < half_num_jobs_node1 + 1; ++i ) {
			auto job = jobgen.register_new_job( 1, i, 1 );
			TS_ASSERT( jobgen._get_job_node( 1, i ) );
			TS_ASSERT_EQUALS( job, jobgen._get_job_node( 1, i ) );
			TS_ASSERT_EQUALS( jobgen.input_source_for_job( 1, i ), 1 );

			jobgen.register_new_job( 1, half_num_jobs_node1 + i, 2 );
			TS_ASSERT_EQUALS( jobgen.input_source_for_job( 1, half_num_jobs_node1 + i ), 2 );
		}

		for ( unsigned int i = 1; i <= num_jobs_node1; ++i ) {
			jobgen.note_job_completed( 1, i, 0 );
		}

		std::list< std::pair< core::Size, core::Size > > container_for_output;
		jobgen.all_job_results_for_node( 1, container_for_output );
		TS_ASSERT_EQUALS( container_for_output.size(), 0 );

		std::list< std::pair< core::Size, core::Size > > container_for_discarded_result_ids;
		jobgen.garbage_collection( 1, true, container_for_discarded_result_ids );
		TS_ASSERT_EQUALS( container_for_discarded_result_ids.size(), 0 );

	}

	void test_job_genealogist2_one_result_for_node1(){
		TR << "beginning test_job_genealogist2_one_result_for_node1()" << std::endl;

		unsigned int const num_job_dag_nodes = 3;
		unsigned int const num_input_sources = 2;
		JobGenealogist_ jobgen( num_job_dag_nodes, num_input_sources );

		unsigned int const num_jobs_node1 = 10;
		unsigned int const half_num_jobs_node1 = num_jobs_node1/2;
		for ( unsigned int i = 1; i < half_num_jobs_node1+1; ++i ) {
			jobgen.register_new_job( 1, i, 1 );
			TS_ASSERT_EQUALS( jobgen.input_source_for_job( 1, i ), 1 );

			jobgen.register_new_job( 1, half_num_jobs_node1 + i, 2 );
			TS_ASSERT_EQUALS( jobgen.input_source_for_job( 1, half_num_jobs_node1 + i ), 2 );
		}

		for ( unsigned int i = 1; i <= num_jobs_node1; ++i ) {
			jobgen.note_job_completed( 1, i, 1 );
		}

		std::list< std::pair< core::Size, core::Size > > container_for_output;
		jobgen.all_job_results_for_node( 1, container_for_output );
		TS_ASSERT_EQUALS( container_for_output.size(), 10 );

		unsigned int count = 0;
		for ( std::pair< core::Size, core::Size > const & ipair : container_for_output ) {
			++count;
			TS_ASSERT_EQUALS( ipair.first, count );
			TS_ASSERT_EQUALS( ipair.second, 1 );
		}

		std::list< std::pair< core::Size, core::Size > > container_for_discarded_result_ids;
		jobgen.garbage_collection( 1, true, container_for_discarded_result_ids );
		TS_ASSERT_EQUALS( container_for_discarded_result_ids.size(), 10 );
		count = 0;
		for ( std::pair< core::Size, core::Size > const & ipair : container_for_output ) {
			++count;
			TS_ASSERT_EQUALS( ipair.first, count );
			TS_ASSERT_EQUALS( ipair.second, 1 );
		}

		std::list< std::pair< core::Size, core::Size > > container_for_output2;
		jobgen.all_job_results_for_node( 1, container_for_output2 );
		TS_ASSERT_EQUALS( container_for_output2.size(), 0 );

	}

	void test_job_genealogist2_simple_test1(){

		/* TREE
		1    2    3    4    5    6    7    8    9    10
		|    |    |    |    |    |    |    |    |    |
		R    R    R    R    R    R    R    R    R    R
		|         |         |         |         |
		11        12        13        14        15
		|         |         |         |         |
		R         R         R         R         R
		|
		16
		|
		R

		*/
		TR << "beginning test_job_genealogist2_simple_test1()" << std::endl;

		unsigned int const num_job_dag_nodes = 3;
		unsigned int const num_input_sources = 2;
		JobGenealogist_ jobgen( num_job_dag_nodes, num_input_sources );

		unsigned int const num_jobs_node1 = 10;
		unsigned int const half_num_jobs_node1 = num_jobs_node1 / 2;
		for ( unsigned int i = 1; i <= num_jobs_node1; ++i ) {
			core::Size input_source_id = ( i <= half_num_jobs_node1 ? 1 : 2 );
			jobgen.register_new_job( 1, i, input_source_id );
			TS_ASSERT_EQUALS( jobgen.input_source_for_job( 1, i ), input_source_id );
		}

		for ( unsigned int i=1; i <= num_jobs_node1; ++i ) {
			jobgen.note_job_completed( 1, i, 1 );
		}

		std::list< std::pair< core::Size, core::Size > > container_for_output;
		jobgen.all_job_results_for_node( 1, container_for_output );
		TS_ASSERT_EQUALS( container_for_output.size(), num_jobs_node1 );
		unsigned int count = 0;
		for ( std::pair< core::Size, core::Size > const & ipair : container_for_output ) {
			++count;
			TS_ASSERT_EQUALS( ipair.first, count );
			TS_ASSERT_EQUALS( ipair.second, 1 );
		}

		/////////
		//NODE 2

		unsigned int const num_jobs_node2 = 5;
		for ( unsigned int i = 1; i <= num_jobs_node2; ++i ) {
			auto const global_job_id = num_jobs_node1 + i;
			auto const parent_global_job_id = 2 * i;
			jobgen.register_new_job( 2, global_job_id, 1, parent_global_job_id, 1 );
			TS_ASSERT_EQUALS( jobgen.input_source_for_job( 2, global_job_id ), jobgen.input_source_for_job( 1, parent_global_job_id ) );
		}

		for ( unsigned int i = 1; i <= num_jobs_node2; ++i ) {
			jobgen.note_job_completed( 2, num_jobs_node1 + i, 1 );
		}

		std::list< std::pair< core::Size, core::Size > > container_for_output2;
		jobgen.all_job_results_for_node( 2, container_for_output2 );
		TS_ASSERT_EQUALS( container_for_output2.size(), num_jobs_node2 );
		count = 0;
		for ( std::pair< core::Size, core::Size > const & ipair : container_for_output2 ) {
			++count;
			TS_ASSERT_EQUALS( ipair.first, num_jobs_node1 + count );
			TS_ASSERT_EQUALS( ipair.second, 1 );
		}

		//////////////////////////////
		//GARBAGE COLLECTION ON NODE 1
		std::list< std::pair< core::Size, core::Size > > container_for_discarded_result_ids;
		jobgen.garbage_collection( 1, true, container_for_discarded_result_ids );
		TS_ASSERT_EQUALS( container_for_discarded_result_ids.size(), num_jobs_node1 - num_jobs_node2 );
		count = 1;
		for ( std::pair< core::Size, core::Size > const & ipair : container_for_discarded_result_ids ) {
			TS_ASSERT_EQUALS( ipair.first, count );
			TS_ASSERT_EQUALS( ipair.second, 1 );
			count += 2;
		}

		std::list< std::pair< core::Size, core::Size > > container_for_output1_again;
		jobgen.all_job_results_for_node( 1, container_for_output1_again );
		TS_ASSERT_EQUALS( container_for_output1_again.size(), num_jobs_node1 - ( num_jobs_node1 - num_jobs_node2 ) );
		count = 2;
		for ( std::pair< core::Size, core::Size > const & ipair : container_for_output1_again ) {
			TS_ASSERT_EQUALS( ipair.first, count );
			TS_ASSERT_EQUALS( ipair.second, 1 );
			count += 2;
		}

		////////
		//NODE 3
		TS_ASSERT(
			jobgen.register_new_job(
			3, //node_id
			num_jobs_node1 + num_jobs_node2 + 1, //global_job_id
			2, //node_id_of_parent
			num_jobs_node1 + 1, //global_job_id_of_parent
			1 //result_id_of_parent
			)
		);
		jobgen.note_job_completed( 3, num_jobs_node1 + num_jobs_node2 + 1, 1 );

		///////////////////////////////////
		//GARBAGE COLLECTION ON NODES 1 & 2
		std::list< std::pair< core::Size, core::Size > > container_for_discarded_result_ids2;
		jobgen.garbage_collection( 2, true, container_for_discarded_result_ids2 );
		jobgen.garbage_collection( 1, true, container_for_discarded_result_ids2 );
		TS_ASSERT_EQUALS( container_for_discarded_result_ids2.size(), 8 );

		std::list< std::pair< core::Size, core::Size > > container_for_output_all;
		jobgen.all_job_results_for_node( 1, container_for_output_all );
		jobgen.all_job_results_for_node( 2, container_for_output_all );
		jobgen.all_job_results_for_node( 3, container_for_output_all );
		TS_ASSERT_EQUALS( container_for_output_all.size(), 3 );

		count = 0;
		for ( std::pair< core::Size, core::Size > const & ipair : container_for_output_all ) {
			switch( ++count ){
			case( 1 ) :
				TS_ASSERT_EQUALS( ipair.first, 2 );
				TS_ASSERT_EQUALS( ipair.second, 1 );
				break;
			case( 2 ) :
				TS_ASSERT_EQUALS( ipair.first, 11 );
				TS_ASSERT_EQUALS( ipair.second, 1 );
				break;
			case( 3 ) :
				TS_ASSERT_EQUALS( ipair.first, 16 );
				TS_ASSERT_EQUALS( ipair.second, 1 );
				break;
			}
		}

		jobgen.print_all_nodes();

	}

	void test_newick_tree(){

		/*
		1     2     3     4     5     6     7     8     9     10
		|     |     |     |     |     |     |     |     |     |
		R     R     R     R     R     R     R     R     R     R
		|     |     |     |     |     |     |     |     |     |
		11    12    13    14    15    16    17    18    19    20
		|     |     |     |     |     |     |     |     |     |
		|\    |\    |\    |\    |\    |\    |\    |\    |\    |\
		| \   | \   | \   | \   | \   | \   | \   | \   | \   | \
		|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
		R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R
		|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
		21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
		|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
		x  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R

		^ x gets deleted at the end
		*/

		TR << "beginning test_newick_tree()" << std::endl;
		unsigned int const num_job_dag_nodes = 3;
		unsigned int const num_input_sources = 2;
		JobGenealogist_ jobgen( num_job_dag_nodes, num_input_sources );

		for ( core::Size i = 1; i <= 10; ++i ) {
			jobgen.register_new_job( 1, i,
				( i < 6 ? 1 : 2 )
			);
			jobgen.note_job_completed( 1, i, 1 );
		}

		for ( core::Size i = 1; i <= 10; ++i ) {
			auto const job_node_ptr1 = jobgen.register_new_job(
				2, 10 + i,
				1, i, 1
			);
			TS_ASSERT( job_node_ptr1 );

			jobgen.note_job_completed( 2, 10 + i, 2 );
			auto const job_node_ptr2 = jobgen._get_job_node( 2, 10 + i );
			TS_ASSERT_EQUALS( job_node_ptr1, job_node_ptr2 );
			TS_ASSERT_EQUALS( job_node_ptr1->children().size(), 2 );
		}

		core::Size count = 0;
		for ( core::Size i = 1; i <= 10; ++i ) {
			core::Size id = 20 + ++count;
			jobgen.register_new_job(
				3, id,
				2, 10 + i, 1
			);

			id = 20 + ++count;
			jobgen.register_new_job(
				3, id,
				2, 10 + i, 2
			);
		}

		for ( core::Size i = 1; i <= 20; ++i ) {
			jobgen.note_job_completed( 3, 20 + i, 1 );
		}

		jobgen.discard_job_result( 3, 21, 1 );//

		std::list< std::pair< core::Size, core::Size > > container_for_output;
		jobgen.garbage_collection( 2, true, container_for_output );
		jobgen.garbage_collection( 1, true, container_for_output );
		TS_ASSERT_EQUALS( container_for_output.size(), 1 );
		TS_ASSERT_EQUALS( container_for_output.begin()->first, 11 );
		TS_ASSERT_EQUALS( container_for_output.begin()->second, 1 );
		/*container_for_output.erase( container_for_output.begin() );
		for ( auto const & p : container_for_output ) {
		TR << "Bad Element: " << p.first << " " << p.second << std::endl;
		}*/

		std::string const tree = jobgen.newick_tree();
		TR << "tree: " << tree << std::endl;
		TS_ASSERT_EQUALS( tree, "((((JR_22_1)JR_11_2)JR_1_1,((JR_23_1)JR_12_1,(JR_24_1)JR_12_2)JR_2_1,((JR_25_1)JR_13_1,(JR_26_1)JR_13_2)JR_3_1,((JR_27_1)JR_14_1,(JR_28_1)JR_14_2)JR_4_1,((JR_29_1)JR_15_1,(JR_30_1)JR_15_2)JR_5_1)input_source_1,(((JR_31_1)JR_16_1,(JR_32_1)JR_16_2)JR_6_1,((JR_33_1)JR_17_1,(JR_34_1)JR_17_2)JR_7_1,((JR_35_1)JR_18_1,(JR_36_1)JR_18_2)JR_8_1,((JR_37_1)JR_19_1,(JR_38_1)JR_19_2)JR_9_1,((JR_39_1)JR_20_1,(JR_40_1)JR_20_2)JR_10_1)input_source_2)all" );
	}

};
