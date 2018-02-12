// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/multistage_rosetta_scripts/MRSJobQueen.cxxtest.hh
/// @brief  test suite for the MRSJobQueen and MRSJob classes
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
//#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/multistage_rosetta_scripts/MRSJobQueen.hh>
#include <protocols/multistage_rosetta_scripts/MRSJob.hh>
#include <protocols/multistage_rosetta_scripts/TagManager.hh>

// Package headers
#include <protocols/jd3/standard/StandardJobQueen.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>

// Movers and Filters
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/filters/BasicFilters.hh>//True Filter

// basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// core headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/options/keys/OptionKey.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/tag/Tag.hh>

using namespace utility::tag;
using namespace protocols::jd3;
using namespace protocols::jd3::standard;
using namespace protocols::multistage_rosetta_scripts;

static basic::Tracer TR("protocols.multistage_rosetta_scripts.MRSJobQueen.cxxtest");

class ReturnInputPoseTenTimesMover : public protocols::relax::FastRelax {

public:
	ReturnInputPoseTenTimesMover( int num_times_to_return ){
		init();
		num_times_to_return_ = num_times_to_return;
	}

	~ReturnInputPoseTenTimesMover(){}

	void apply( core::pose::Pose & pose ) override {
		//TR << "Running ReturnInputPoseTenTimesMover::apply()" <<std::endl;
		pose_ = pose.clone();
		return_count_ = 1;
		apply_was_called_ = true;
	}

	bool apply_was_called() const { return apply_was_called_; }

	core::pose::PoseOP get_additional_output() override {
		if ( ++return_count_ <= num_times_to_return_ ) {
			return pose_;
		} else {
			return 0;
		}
	}

private:
	void init(){
		return_count_ = 0;
		apply_was_called_ = false;
	}

	int return_count_;
	core::pose::PoseOP pose_;
	bool apply_was_called_;

	int num_times_to_return_;
};

class MRSJobQueenTests : public CxxTest::TestSuite
{
public:

	MRSJobQueenTests() {}

	void setUp() override {
		core_init();

		job_def_filename_ = "protocols/multistage_rosetta_scripts/job_def_file.xml" ;
		pose1_filename_ = "protocols/multistage_rosetta_scripts/3U3B_A.pdb" ;
		pose2_filename_ = "protocols/multistage_rosetta_scripts/3U3B_B.pdb" ;

		result_for_every_job_.reset( new core::pose::Pose );
		core::import_pose::pose_from_file( *result_for_every_job_, pose1_filename_, core::import_pose::PDB_file);

		result_vec_.resize( 1 );
		result_vec_[ 1 ].reset( new protocols::jd3::standard::PoseJobResult( result_for_every_job_ ) );
	}

	void test_MRSJobQueen_with_good_input() {
		TS_ASSERT( result_for_every_job_ );
		//basic::options::option[ basic::options::OptionKeys::out::keep_ancestor_poses ]( true );
		basic::options::option[ basic::options::OptionKeys::in::file::job_definition_file ]( job_def_filename_ );
		MRSJobQueen zero_queen;
		MRSJobQueen worker_queen;

		JobDigraphOP job_dag = zero_queen.initial_job_dag();
		check_job_dag( job_dag );

		inner_test_node1( zero_queen, worker_queen );
		inner_test_node2( zero_queen, worker_queen );
		inner_test_node3( zero_queen, worker_queen );

		//TODO The code has changed enough since this test was written that this is no longer valid. Jack needs to do the work to fix this
		//inner_test_node4( zero_queen, worker_queen );

		inner_test_ReturnInputPoseTenTimesMover();
	}

	static void check_job_dag( JobDigraphOP dag ){
		//Nodes
		TS_ASSERT_EQUALS( dag->num_nodes(), 4 );

		//Edges
		TS_ASSERT_EQUALS( dag->num_edges(), 3 );
		TS_ASSERT( dag->find_job_edge( 1, 2 ) );
		TS_ASSERT( dag->find_job_edge( 2, 3 ) );
		TS_ASSERT( dag->find_job_edge( 3, 4 ) );
	}

	void inner_test_node1( MRSJobQueen & zero_queen, MRSJobQueen & worker_queen ) {
		//should be a total of 10 jobs. Let's collect them 4 at a time
		std::list< LarvalJobOP > ljobs_1_through_4 = zero_queen.determine_job_list( 1, 4 );
		TS_ASSERT_EQUALS( ljobs_1_through_4.size(), 4 );

		std::list< LarvalJobOP > ljobs_5_through_8 = zero_queen.determine_job_list( 1, 4 );
		TS_ASSERT_EQUALS( ljobs_5_through_8.size(), 4 );

		std::list< LarvalJobOP > ljobs_9_through_10 = zero_queen.determine_job_list( 1, 4 );
		TS_ASSERT_EQUALS( ljobs_9_through_10.size(), 2 );

		std::list< LarvalJobOP > hopefully_no_ljobs = zero_queen.determine_job_list( 1, 4 );
		TS_ASSERT_EQUALS( hopefully_no_ljobs.size(), 0 );

		utility::vector1< LarvalJobOP > all_ljobs;
		all_ljobs.reserve( 10 );
		for ( LarvalJobOP ljob : ljobs_1_through_4 ) all_ljobs.push_back( ljob );
		for ( LarvalJobOP ljob : ljobs_5_through_8 ) all_ljobs.push_back( ljob );
		for ( LarvalJobOP ljob : ljobs_9_through_10 ) all_ljobs.push_back( ljob );
		TS_ASSERT_EQUALS( all_ljobs.size(), 10 );

		for ( int i = 2; i < 6; ++i ) {
			TS_ASSERT_EQUALS( all_ljobs[ 1 ]->inner_job(), all_ljobs[ i ]->inner_job() );
			TS_ASSERT_EQUALS( all_ljobs[ 6 ]->inner_job(), all_ljobs[ 5 + i ]->inner_job() );

			TS_ASSERT_EQUALS( static_cast< StandardInnerLarvalJob const & > ( * all_ljobs[ i ]->inner_job() ).prelim_job_node(), 1 );
			TS_ASSERT_EQUALS( static_cast< StandardInnerLarvalJob const & > ( * all_ljobs[ 5 + i ]->inner_job() ).prelim_job_node(), 2 );
		}

		utility::vector1< JobResultCOP > dummy;

		for ( LarvalJobOP ljob : all_ljobs ) {
			JobOP job_op = worker_queen.mature_larval_job( ljob, dummy );
			//JobOP job_op = worker_queen.complete_larval_job_maturation( ljob, 0, dummy );
			MRSJob & job = static_cast< MRSJob & > ( * job_op );
			std::list< mover_or_filter > & protocols = job.protocols();

			TS_ASSERT_EQUALS( protocols.size(), 2 );

			TS_ASSERT( protocols.front().mover );
			TS_ASSERT( ! protocols.front().filter );
			TS_ASSERT_EQUALS( protocols.front().mover->get_name(), protocols::relax::FastRelax::mover_name() );

			TS_ASSERT( ! protocols.back().mover );
			TS_ASSERT( protocols.back().filter );
			TS_ASSERT_EQUALS( protocols.back().filter->name(), protocols::simple_filters::ScoreTypeFilter::class_name() );

			TS_ASSERT_EQUALS( zero_queen.stage_for_global_job_id( ljob->job_index() ), 1 );
			TS_ASSERT_EQUALS( worker_queen.stage_for_global_job_id( ljob->job_index() ) + ljob->job_index(), 1 + ljob->job_index() );
		}

		//return results
		std::list< core::Real > job_result_scores;
		job_result_scores.push_back( -4  );//1,1
		job_result_scores.push_back( -22 );//2,1
		job_result_scores.push_back( -80 );//2,2 KEEP
		job_result_scores.push_back( -49 );//3,1
		job_result_scores.push_back( -1  );//4,1
		job_result_scores.push_back( -50 );//4,2 KEEP
		job_result_scores.push_back( -4  );//5,1

		job_result_scores.push_back( -16 );//6,1 KEEP
		job_result_scores.push_back( -2  );//6,2
		job_result_scores.push_back( -6  );//7,1 KEEP
		job_result_scores.push_back( -4  );//8,1
		job_result_scores.push_back( -1  );//8,2
		job_result_scores.push_back( -1  );//9,1
		job_result_scores.push_back( -4  );//10,1
		job_result_scores.push_back( -4  );//10,2


		std::list< core::Real >::iterator it = job_result_scores.begin();

		for ( core::Size global_job_id = 1; global_job_id <= 10; ++global_job_id ) {
			//2 results for even numbers, 1 results for odd numbers
			core::Size const num_results = 1 + ( global_job_id % 2 == 0 ? 1 : 0 );
			zero_queen.note_job_completed( global_job_id, jd3_job_status_success, num_results, true );

			for ( core::Size job_result_id = 1; job_result_id <= num_results; ++job_result_id ) {
				JobSummaryOP summary( new protocols::jd3::standard::EnergyJobSummary( *it ) );
				++it;
				zero_queen.completed_job_summary( global_job_id, job_result_id, summary );
			}
		}

		std::list< std::pair< core::Size, core::Size > > jobs_that_should_be_discarded = zero_queen.job_results_that_should_be_discarded();
		TS_ASSERT_EQUALS( jobs_that_should_be_discarded.size(), job_result_scores.size() - 4 );

		//Check for absence of good results
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 2, 2 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 4, 2 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 6, 1 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 7, 1 ) ), jobs_that_should_be_discarded.end() );

		//Check for presence of bad results
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 1, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 2, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 3, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 4, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 5, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 6, 2 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 8, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 8, 2 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 9, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 10, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 10, 2 ) ) != jobs_that_should_be_discarded.end() );

	}

	void inner_test_node2( MRSJobQueen & zero_queen, MRSJobQueen & worker_queen ) {
		//make sure no jobs are given for node 1
		TS_ASSERT_EQUALS( zero_queen.determine_job_list( 1, 4 ).size(), 0 );

		//should be a total of 8 jobs. Let's collect them 4 at a time
		std::list< LarvalJobOP > ljobs_1_through_4 = zero_queen.determine_job_list( 2, 4 );
		TS_ASSERT_EQUALS( ljobs_1_through_4.size(), 4 );

		std::list< LarvalJobOP > ljobs_5_through_8 = zero_queen.determine_job_list( 2, 4 );
		TS_ASSERT_EQUALS( ljobs_5_through_8.size(), 4 );

		std::list< LarvalJobOP > hopefully_no_ljobs = zero_queen.determine_job_list( 2, 4 );
		TS_ASSERT_EQUALS( hopefully_no_ljobs.size(), 0 );

		utility::vector1< LarvalJobOP > all_ljobs;
		all_ljobs.reserve( 8 );
		for ( LarvalJobOP ljob : ljobs_1_through_4 ) all_ljobs.push_back( ljob );
		for ( LarvalJobOP ljob : ljobs_5_through_8 ) all_ljobs.push_back( ljob );
		TS_ASSERT_EQUALS( all_ljobs.size(), 8 );

		for ( int i = 1; i <= 4; ++i ) {
			TS_ASSERT_EQUALS( static_cast< StandardInnerLarvalJob const & > ( * all_ljobs[ i ]->inner_job() ).prelim_job_node(), 1 );
			TS_ASSERT_EQUALS( static_cast< StandardInnerLarvalJob const & > ( * all_ljobs[ 4 + i ]->inner_job() ).prelim_job_node(), 2 );
		}

		core::Size local_job_id = 0;
		for ( LarvalJobOP ljob : all_ljobs ) {
			++local_job_id;
			JobOP job_op = worker_queen.mature_larval_job( ljob, result_vec_ );
			//JobOP job_op = worker_queen.complete_larval_job_maturation( ljob, 0, result_vec_ );

			core::Size const global_job_id = local_job_id + 10;
			TS_ASSERT_EQUALS( global_job_id, ljob->job_index() );

			MRSJob & job = static_cast< MRSJob & > ( * job_op );
			std::list< mover_or_filter > & protocols = job.protocols();

			TS_ASSERT_EQUALS( protocols.size(), 3 );

			TS_ASSERT( protocols.front().mover );
			TS_ASSERT( ! protocols.front().filter );
			TS_ASSERT_EQUALS( protocols.front().mover->get_name(), protocols::relax::FastRelax::mover_name() );

			std::list< mover_or_filter >::const_iterator middle = std::next( protocols.begin() );
			TS_ASSERT( middle->mover );
			TS_ASSERT( ! middle->filter );
			TS_ASSERT_EQUALS( middle->mover->get_name(), protocols::idealize::IdealizeMover::mover_name() );

			TS_ASSERT( ! protocols.back().mover );
			TS_ASSERT( protocols.back().filter );
			TS_ASSERT_EQUALS( protocols.back().filter->name(), protocols::simple_filters::ScoreTypeFilter::class_name() );

			TS_ASSERT_EQUALS( zero_queen.stage_for_global_job_id( ljob->job_index() ), 2 );
			TS_ASSERT_EQUALS( worker_queen.stage_for_global_job_id( ljob->job_index() ) + ljob->job_index(), 2 + ljob->job_index() );
		}

		//return results
		//Returning 10 results that will be segregated into 2 segregation, 5 each.
		std::list< core::Real > job_result_scores;
		job_result_scores.push_back( -4  );//11,1
		job_result_scores.push_back( -80 );//12,1 KEEP
		job_result_scores.push_back( -70 );//12,2 KEEP

		job_result_scores.push_back( -9  );//13,1
		job_result_scores.push_back( -1  );//14,1
		job_result_scores.push_back( -50 );//14,2 KEEP

		job_result_scores.push_back( -4  );//15,1
		job_result_scores.push_back( -6  );//16,1
		job_result_scores.push_back( -2  );//16,2

		job_result_scores.push_back( -40 );//17,1 KEEP
		job_result_scores.push_back( -4  );//18,1
		job_result_scores.push_back( -1  );//18,2

		std::list< core::Real >::iterator it = job_result_scores.begin();
		for ( core::Size global_job_id = 11; global_job_id <= 18; ++global_job_id ) {
			//2 results for even numbers, 1 results for odd numbers
			core::Size const num_results = 1 + ( global_job_id % 2 == 0 ? 1 : 0 );
			zero_queen.note_job_completed( global_job_id, jd3_job_status_success, num_results, true );

			for ( core::Size job_result_id = 1; job_result_id <= num_results; ++job_result_id ) {
				JobSummaryOP summary( new protocols::jd3::standard::EnergyJobSummary( *it ) );
				++it;
				zero_queen.completed_job_summary( global_job_id, job_result_id, summary );
			}
		}

		std::list< std::pair< core::Size, core::Size > > jobs_that_should_be_discarded = zero_queen.job_results_that_should_be_discarded();
		TS_ASSERT_EQUALS( jobs_that_should_be_discarded.size(), job_result_scores.size() - 4 + 1 );//1 node 1 job result can be discarded
		//We expect (6,1) to also be discarded here

		//Check for absence of good results
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 12, 1 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 12, 2 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 14, 2 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 17, 1 ) ), jobs_that_should_be_discarded.end() );

		//Check for presence of bad results
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 11, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 13, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 14, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 15, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 16, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 16, 2 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 18, 1 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 18, 2 ) ) != jobs_that_should_be_discarded.end() );

		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 6, 1 ) ) != jobs_that_should_be_discarded.end() );
	}

	void inner_test_node3( MRSJobQueen & zero_queen, MRSJobQueen & worker_queen ) {

		TS_ASSERT_EQUALS( zero_queen.determine_job_list( 1, 1 ).size(), 0 );
		TS_ASSERT_EQUALS( zero_queen.determine_job_list( 2, 1 ).size(), 0 );

		//should be a total of 4 jobs. Let's collect them 4 at a time
		std::list< LarvalJobOP > ljobs_1_through_4 = zero_queen.determine_job_list( 3, 4 );
		TS_ASSERT_EQUALS( ljobs_1_through_4.size(), 4 );

		std::list< LarvalJobOP > hopefully_no_ljobs = zero_queen.determine_job_list( 3, 4 );
		TS_ASSERT_EQUALS( hopefully_no_ljobs.size(), 0 );

		utility::vector1< LarvalJobOP > all_ljobs;
		all_ljobs.reserve( 4 );
		for ( LarvalJobOP ljob : ljobs_1_through_4 ) all_ljobs.push_back( ljob );
		TS_ASSERT_EQUALS( all_ljobs.size(), 4 );

		core::Size local_job_id = 0;
		for ( LarvalJobOP ljob : all_ljobs ) {
			++local_job_id;
			JobOP job_op = worker_queen.mature_larval_job( ljob, result_vec_ );

			core::Size const global_job_id = local_job_id + 18;
			TS_ASSERT_EQUALS( global_job_id, ljob->job_index() );

			MRSJob & job = static_cast< MRSJob & > ( * job_op );
			std::list< mover_or_filter > & protocols = job.protocols();

			TS_ASSERT_EQUALS( protocols.size(), 3 );

			TS_ASSERT( protocols.front().mover );
			TS_ASSERT( ! protocols.front().filter );
			TS_ASSERT_EQUALS( protocols.front().mover->get_name(), protocols::relax::FastRelax::mover_name() );

			std::list< mover_or_filter >::const_iterator middle = std::next( protocols.begin() );
			TS_ASSERT( ! middle->mover );
			TS_ASSERT( middle->filter );
			TS_ASSERT_EQUALS( middle->filter->name(), protocols::simple_filters::InterfaceSasaFilter::class_name() );

			TS_ASSERT( ! protocols.back().mover );
			TS_ASSERT( protocols.back().filter );
			TS_ASSERT_EQUALS( protocols.back().filter->name(), protocols::simple_filters::ScoreTypeFilter::class_name() );

			TS_ASSERT_EQUALS( zero_queen.stage_for_global_job_id( ljob->job_index() ), 3 );
			TS_ASSERT_EQUALS( worker_queen.stage_for_global_job_id( ljob->job_index() ) + ljob->job_index(), 3 + ljob->job_index() );
		}

		//return results
		//Returning 10 results that will be segregated into 2 segregation, 5 each.
		std::list< core::Real > job_result_scores;
		job_result_scores.push_back( -95 );//19,1
		job_result_scores.push_back( -90 );//19,2
		job_result_scores.push_back( -4  );//19,3
		job_result_scores.push_back( -4  );//19,4

		//20 has no results

		job_result_scores.push_back( -85 );//21,1
		job_result_scores.push_back( -4  );//21,2
		job_result_scores.push_back( -4  );//21,3

		job_result_scores.push_back( -4  );//22,1
		job_result_scores.push_back( -4  );//22,2
		job_result_scores.push_back( -4  );//22,3
		job_result_scores.push_back( -99 );//22,4
		job_result_scores.push_back( -4  );//22,5

		zero_queen.note_job_completed( 20, jd3_job_status_failed_do_not_retry, 99, true );
		//JobSummaryOP summary20( new protocols::jd3::standard::EnergyJobSummary( 9 ) );

		std::list< core::Real >::iterator it = job_result_scores.begin();
		for ( core::Size global_job_id = 19; global_job_id <= 22; ++global_job_id ) {
			if ( global_job_id == 20 ) continue;

			//2 results for even numbers, 1 results for odd numbers
			core::Size num_results = 0;
			switch( global_job_id ){
			case( 19 ) :
				num_results = 4;
				break;
			case( 21 ) :
				num_results = 3;
				break;
			case( 22 ) : //TODO access node managers here?
				num_results = 5;
				break;
			}

			zero_queen.note_job_completed( global_job_id, jd3_job_status_success, num_results, true );

			for ( core::Size job_result_id = 1; job_result_id <= num_results; ++job_result_id ) {
				JobSummaryOP summary( new protocols::jd3::standard::EnergyJobSummary( *it ) );
				++it;
				zero_queen.completed_job_summary( global_job_id, job_result_id, summary );
			}
		}

		std::list< std::pair< core::Size, core::Size > > jobs_that_should_be_discarded = zero_queen.job_results_that_should_be_discarded();
		//for ( std::pair< core::Size, core::Size > waste : jobs_that_should_be_discarded ) {
		//std::cout << "mmm " << waste.first << " " << waste.second << std::endl;
		//}
		TS_ASSERT_EQUALS( jobs_that_should_be_discarded.size(), 12 );//1 node 1 job result can be discarded + 2 node 2 + 9 node 3

		//TODO the treatment of these results have since changed since this test has been designed. Jack needs to do the work to figure out what numbers should be here now.
		/*
		//Check for absence of good results
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 19, 1 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 19, 2 ) ), jobs_that_should_be_discarded.end() );
		TS_ASSERT_EQUALS( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 20, 1 ) ), jobs_that_should_be_discarded.end() );

		//Check for presence of bad results
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 19, 3 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 19, 4 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 21, 2 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 21, 3 ) ) != jobs_that_should_be_discarded.end() );

		//node 2
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 12, 2 ) ) != jobs_that_should_be_discarded.end() );
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 17, 1 ) ) != jobs_that_should_be_discarded.end() );

		//node 1
		TS_ASSERT( std::find( jobs_that_should_be_discarded.begin(), jobs_that_should_be_discarded.end(), std::pair< core::Size, core::Size >( 7, 1 ) ) != jobs_that_should_be_discarded.end() );
		*/
	}


	void inner_test_node4( MRSJobQueen & zero_queen, MRSJobQueen & worker_queen ) {

		TS_ASSERT_EQUALS( zero_queen.determine_job_list( 1, 1 ).size(), 0 );
		TS_ASSERT_EQUALS( zero_queen.determine_job_list( 2, 1 ).size(), 0 );
		TS_ASSERT_EQUALS( zero_queen.determine_job_list( 3, 1 ).size(), 0 );

		std::list< LarvalJobOP > ljobs_1_through_3 = zero_queen.determine_job_list( 4, 4 );
		TS_ASSERT_EQUALS( ljobs_1_through_3.size(), 3 );

		std::list< LarvalJobOP > hopefully_no_ljobs = zero_queen.determine_job_list( 4, 4 );
		TS_ASSERT_EQUALS( hopefully_no_ljobs.size(), 0 );

		core::Size local_job_id = 0;
		for ( LarvalJobOP ljob : ljobs_1_through_3 ) {
			++local_job_id;
			JobOP job_op = worker_queen.mature_larval_job( ljob, result_vec_ );
			//JobOP job_op = worker_queen.complete_larval_job_maturation( ljob, 0, result_vec_ );
			MRSJob & job = static_cast< MRSJob & > ( * job_op );
			TS_ASSERT( job.pose() );

			std::list< mover_or_filter > & protocols = job.protocols();

			TS_ASSERT_EQUALS( job.max_num_results(), 3 );
			TS_ASSERT_EQUALS( protocols.size(), 9 );

			//FIRST
			for ( mover_or_filter & mof : protocols ) {
				if ( mof.mover ) {
					TS_ASSERT( ! mof.filter );
					mof.mover.reset( new ReturnInputPoseTenTimesMover( 50 ) );
				} else {
					TS_ASSERT( mof.filter );
					mof.filter.reset( new protocols::filters::TrueFilter );
				}
			}

			utility::pointer::static_pointer_cast< MRSJob >( job_op )->set_pose( result_for_every_job_->clone() );
			TS_ASSERT( utility::pointer::static_pointer_cast< MRSJob >( job_op )->pose() );
			CompletedJobOutput cjo = job_op->run();
			TS_ASSERT_EQUALS( cjo.job_results.size(), 3 );
			TS_ASSERT_EQUALS( cjo.status, jd3_job_status_success );

			//SECOND
			for ( mover_or_filter & mof : protocols ) {
				if ( mof.mover ) {
					mof.mover.reset( new ReturnInputPoseTenTimesMover( 2 ) );
				} else {
					mof.filter.reset( new protocols::filters::TrueFilter );
				}
			}

			utility::pointer::static_pointer_cast< MRSJob >( job_op )->set_pose( result_for_every_job_->clone() );
			TS_ASSERT( utility::pointer::static_pointer_cast< MRSJob >( job_op )->pose() );
			cjo = job_op->run();
			TS_ASSERT_EQUALS( cjo.job_results.size(), 2 );
			TS_ASSERT_EQUALS( cjo.status, jd3_job_status_success );

			//THIRD
			for ( mover_or_filter & mof : protocols ) {
				if ( mof.mover ) {
					mof.mover.reset( new ReturnInputPoseTenTimesMover( 1 ) );
				} else {
					mof.filter.reset( new protocols::filters::TrueFilter );
				}
			}

			utility::pointer::static_pointer_cast< MRSJob >( job_op )->set_pose( result_for_every_job_->clone() );
			TS_ASSERT( utility::pointer::static_pointer_cast< MRSJob >( job_op )->pose() );
			cjo = job_op->run();
			TS_ASSERT_EQUALS( cjo.job_results.size(), 1 );
			TS_ASSERT_EQUALS( cjo.status, jd3_job_status_success );

			//FOURTH
			for ( mover_or_filter & mof : protocols ) {
				if ( mof.mover ) {
					mof.mover.reset( new ReturnInputPoseTenTimesMover( 50 ) );
				} else {
					mof.filter.reset( new protocols::filters::FalseFilter );
				}
			}

			utility::pointer::static_pointer_cast< MRSJob >( job_op )->set_pose( result_for_every_job_->clone() );
			TS_ASSERT( utility::pointer::static_pointer_cast< MRSJob >( job_op )->pose() );
			cjo = job_op->run();
			TS_ASSERT_EQUALS( cjo.job_results.size(), 0 );
			TS_ASSERT_EQUALS( cjo.status, jd3_job_status_failed_do_not_retry );
		}

	}

	void inner_test_ReturnInputPoseTenTimesMover(){
		utility::pointer::shared_ptr< ReturnInputPoseTenTimesMover > mover( new ReturnInputPoseTenTimesMover( 11 ) );
		TS_ASSERT( result_for_every_job_ );
		TS_ASSERT( ! mover->apply_was_called() );
		mover->apply( * result_for_every_job_ );
		TS_ASSERT( mover->apply_was_called() );
		for ( short int i = 1; i < 11; ++i ) {
			TS_ASSERT( mover->get_additional_output() );
		}
		for ( short int i = 1; i < 100; ++i ) {
			TS_ASSERT( ! mover->get_additional_output() );
		}
		TS_ASSERT( mover->apply_was_called() );
	}

	void test_not_enough_job_results_for_next_node(){
		TS_ASSERT( result_for_every_job_ );
		//basic::options::option[ basic::options::OptionKeys::parser::protocol ]( "protocols/multistage_rosetta_scripts/good_script1.xml" );
		basic::options::option[ basic::options::OptionKeys::in::file::job_definition_file ]( job_def_filename_ );
		MRSJobQueen zero_queen;
		//MRSJobQueen worker_queen;

		JobDigraphOP job_dag = zero_queen.initial_job_dag();
		check_job_dag( job_dag );

		std::list< LarvalJobOP > all_ljobs = zero_queen.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( all_ljobs.size(), 10 );

		//Node 2 needs 4 job results from node 1, only supply 2
		zero_queen.note_job_completed( 1, jd3_job_status_success, 2, true );
		JobSummaryOP summary1( new protocols::jd3::standard::EnergyJobSummary( 1 ) );
		JobSummaryOP summary2( new protocols::jd3::standard::EnergyJobSummary( 2 ) );
		zero_queen.completed_job_summary( 1, 1, summary1 );
		zero_queen.completed_job_summary( 1, 2, summary2 );

		for ( int i=2; i<=10; ++i ) {
			zero_queen.note_job_completed( i, jd3_job_status_failed_do_not_retry, 0, true );
		}

		std::list< std::pair< core::Size, core::Size > > jobs_that_should_be_discarded = zero_queen.job_results_that_should_be_discarded();
		TS_ASSERT_EQUALS( jobs_that_should_be_discarded.size(), 0 );


		//Check Node 2
		all_ljobs = zero_queen.determine_job_list( 2, 100 );
		TS_ASSERT_EQUALS( all_ljobs.size(), 4 );

		//Give enough node 2 results. Node 3 expects 4
		for ( int i=11; i<=14; ++i ) {
			zero_queen.note_job_completed( i, jd3_job_status_success, 2, true );
			JobSummaryOP summary1( new protocols::jd3::standard::EnergyJobSummary( 1 ) );
			JobSummaryOP summary2( new protocols::jd3::standard::EnergyJobSummary( 2 ) );
			zero_queen.completed_job_summary( i, 1, summary1 );
			zero_queen.completed_job_summary( i, 2, summary2 );
		}

		//Check Node 3
		all_ljobs = zero_queen.determine_job_list( 3, 100 );
		TS_ASSERT_EQUALS( all_ljobs.size(), 4 );
	}

	//NEGAVTIVE TESTS
	void test_tag_for_nonmover_nonfilter(){
		basic::options::option[ basic::options::OptionKeys::in::file::job_definition_file ]( "protocols/multistage_rosetta_scripts/job_def_file.bad1.xml" );
		bool exception_thrown = false;

		try{
			MRSJobQueen queen;
			JobDigraphOP job_dag = queen.initial_job_dag();
		} catch( ... ){
			exception_thrown = true;
		}

		TS_ASSERT( exception_thrown );
	}

	void test_typo_in_score_function(){
		basic::options::option[ basic::options::OptionKeys::in::file::job_definition_file ]( "protocols/multistage_rosetta_scripts/job_def_file.bad2.xml" );
		bool exception_thrown = false;

		try{
			MRSJobQueen queen;
			JobDigraphOP job_dag = queen.initial_job_dag();
			queen.determine_validity_of_stage_tags();

			//std::vector< TagListOP > const & tags = queen.tag_manager().tag_list_for_input_pose_id();
			//for( unsigned int i = 1; i < tags.size(); ++i ){
			// TS_ASSERT( tags[ i ] != tags[ 0 ] );
			//}
		} catch( ... ){
			exception_thrown = true;
		}

		TS_ASSERT( exception_thrown );
	}

	void test_incorrect_attribute_value(){
		basic::options::option[ basic::options::OptionKeys::in::file::job_definition_file ]( "protocols/multistage_rosetta_scripts/job_def_file.bad3.xml" );
		bool exception_thrown = false;

		try{
			MRSJobQueen queen;
			JobDigraphOP job_dag = queen.initial_job_dag();
		} catch( ... ){
			exception_thrown = true;
		}

		TS_ASSERT( exception_thrown );
	}

	//There is nothing wrong with this code, except that it is testing a feature that does not exist yet.
	//Please uncomment this when JD3 supports PDBlists in the job definition file
	//-Jack Maguire, Feb 2018
	/*
	void _test_pdblist(){
	basic::options::option[ basic::options::OptionKeys::in::file::job_definition_file ]( "protocols/multistage_rosetta_scripts/job_def_file.pdblist.xml" );

	try{
	MRSJobQueen queen;
	JobDigraphOP job_dag = queen.initial_job_dag();
	TS_ASSERT_EQUALS( queen.num_input_structs(), 3 );

	TagManager const & tag_manager = queen.tag_manager();
	std::vector< TagListOP > const & tag_lists = tag_manager.tag_list_for_input_pose_id();
	TS_ASSERT_EQUALS( tag_lists[ 0 ], tag_lists[ 1 ] );
	TS_ASSERT( tag_lists[ 1 ] != tag_lists[ 2 ] );

	{
	bool sfxn_exists = false;
	for( auto const & tag : * tag_lists[ 2 ] ){
	if( tag->getName() == "SCOREFXNS" ){
	for( auto const & subtag : tag->getTags() ){
	if( subtag->getName() == "ScoreFunction" ){
	if( subtag->getOption< std::string >( "name" ) == "individual_sfxn" ){
	sfxn_exists = true;
	TS_ASSERT_EQUALS( subtag->getOption< std::string >( "weights" ), "beta_nov15_cst.wts" );
	}
	}
	}
	}
	}
	TS_ASSERT( sfxn_exists );
	}
	} catch( utility::excn::Exception const & e ){
	TR << "caught exception: " << std::endl;
	TR << e.msg() << std::endl;
	TS_ASSERT( false );
	}

	}*/


	//GDB stuff
	static void print_list( std::list< std::pair< core::Size, core::Size > > & l );

	static void gdb_print( utility::vector1< JGJobNode * > const & );

private:
	std::string job_def_filename_; // "protocols/multistage_rosetta_scripts/job_def_file.xml";
	std::string pose1_filename_; // "protocols/multistage_rosetta_scripts/3U3B_A.pdb";
	std::string pose2_filename_; // "protocols/multistage_rosetta_scripts/3U3B_B.pdb";

	core::pose::PoseOP result_for_every_job_;
	utility::vector1< JobResultCOP > result_vec_;

	protocols::jd3::standard::PoseJobResultOP pose_result_;
	protocols::jd3::standard::EnergyJobSummaryOP energy_job_summary_;
};

void MRSJobQueenTests::print_list( std::list< std::pair< core::Size, core::Size > > & l ){
	int j=0;
	for ( std::pair< core::Size, core::Size > & i : l ) {
		TR << (++j) << " " << i.first << " " << i.second << std::endl;
	}
}

void MRSJobQueenTests::gdb_print( utility::vector1< JGJobNode * > const & v ){
	for ( unsigned int i = 1; i <= v.size(); ++i ) {
		TR << i << " " << v[i]->global_job_id() << " " << v[i]->node() << std::endl;
	}
}
