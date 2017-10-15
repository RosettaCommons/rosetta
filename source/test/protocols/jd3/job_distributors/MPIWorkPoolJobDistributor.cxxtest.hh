// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/jd3/MPIWorkPoolJobDistributor.cxxtest.hh
/// @brief  test suite for protocols::jd3::job_distributors/MPIWorkPoolJobDistributor using the SimulateMPI utilities
/// @author Andy Watkins (amw579@nyu.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/mpi_funcs.hh>

// Unit headers
#include <protocols/jd3/job_distributors/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/JobSummary.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>

/// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>
//#include <core/pack/scmin/SidechainStateAssignment.hh>
#include <protocols/moves/NullMover.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>

// Utility headers
#include <utility/SimulateMPI.hh>
#include <utility/mpi_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <map>
#include <sstream>
#include <utility>

#ifdef SERIALIZATION
// Utility headers
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

#endif

using namespace utility;
using namespace core;
using namespace protocols::jd3;
using namespace protocols::jd3::job_distributors;
using namespace protocols::jd3::standard;

// shorthand functions, sp and sp1
std::pair< Size, Size > sp( Size first_in, Size second_in ) { return std::make_pair( first_in, second_in ); }
std::pair< Size, Size > sp1( Size first_in ) { return std::make_pair( first_in, Size(1) ); }

class PoolMnPJob : public MoverAndPoseJob
{
public:

	CompletedJobOutput
	run() override {
		CompletedJobOutput output = MoverAndPoseJob::run();
		EnergyJobSummaryOP sum = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( output.job_results[1].first );
		sum->energy( 1234 );
		return output;
	}
};

typedef utility::pointer::shared_ptr< PoolMnPJob > PoolMnPJobOP;

// Dummy JobQueen //
class PoolJQ1 : public protocols::jd3::standard::StandardJobQueen
{
public:
	PoolJQ1() : njobs_( 1 ), job_list_determined_( false ), n_jobs_created_so_far_( 0 ) {}
	~PoolJQ1() {}

	// create a job graph with a single node; skip the SJQ's implementation of this function
	JobDigraphOP initial_job_dag() override {
		JobDigraphOP job_graph( new JobDigraph( 1 ) );
		return job_graph;
	}

	LarvalJobs determine_job_list( Size job_node_index, Size max_njobs ) override
	{
		LarvalJobs jobs;
		if ( job_list_determined_ ) return jobs;

		StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( njobs_, job_node_index ) );
		jobs = expand_job_list( inner_job, max_njobs );
		n_jobs_created_so_far_ += jobs.size();
		if ( n_jobs_created_so_far_ >= njobs_ ) {
			job_list_determined_ = true;
		}
		return jobs;
	}

	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP job,
		utility::options::OptionCollectionCOP,
		utility::vector1< JobResultCOP > const &
	) override
	{
		using namespace protocols::jd3;
		PoolMnPJobOP mature_job( new PoolMnPJob );

		core::pose::PoseOP pose = create_trpcage_ideal_poseop();
		mature_job->pose( pose );
		mature_job->mover( protocols::moves::NullMoverOP( new protocols::moves::NullMover ) );

		jobs_matured_.push_back( job->job_index() );

		return mature_job;
	}

	void note_job_completed( core::Size job_id, JobStatus status, core::Size nresults  ) override { status_[ job_id ] = status; StandardJobQueen::note_job_completed( job_id, status, nresults ); }
	void completed_job_summary( core::Size job_id, core::Size result_index, JobSummaryOP summary ) override { summaries_[ sp(job_id,result_index) ] = summary; }
	void completed_job_result( LarvalJobCOP job, core::Size result_index, JobResultOP job_result ) override { results_[ sp(job->job_index(),result_index) ] = job_result; }

	bool has_job_completed( protocols::jd3::LarvalJobCOP ) override { return false; }
	void mark_job_as_having_begun( protocols::jd3::LarvalJobCOP /*job*/ ) override {/*TEMP*/}

public:
	core::Size njobs_;

	bool job_list_determined_;
	core::Size n_jobs_created_so_far_;

	mutable utility::vector1< core::Size > jobs_matured_;
	std::map< Size,        JobStatus    > status_;
	std::map< JobResultID, JobSummaryOP > summaries_;
	std::map< JobResultID, JobResultOP  > results_;

};

class PoolMnPJob2 : public MoverAndPoseJob
{
public:
	PoolMnPJob2() : throw_( false ), throw_bad_inputs_( false ), status_( jd3_job_status_success ) {}

	CompletedJobOutput
	run() override {
		CompletedJobOutput output = MoverAndPoseJob::run();
		if ( throw_ ) {
			throw utility::excn::EXCN_Msg_Exception( "PoolMnPJob2 exception" );
		} else if ( throw_bad_inputs_ ) {
			throw utility::excn::EXCN_BadInput( "PoolMnPJob2 bad input exception" );
		} else {
			output.status = status_;
		}
		return output;
	}

	bool throw_;
	bool throw_bad_inputs_;
	JobStatus status_;

};

typedef utility::pointer::shared_ptr< PoolMnPJob2 > PoolMnPJob2OP;

class PoolJQ2 : public PoolJQ1
{
public:
	PoolJQ2() : throw_( false ), throw_bad_inputs_( false ), status_( jd3_job_status_success ) {}

	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP job,
		utility::options::OptionCollectionCOP,
		utility::vector1< JobResultCOP > const &
	) override
	{
		using namespace protocols::jd3;
		PoolMnPJob2OP mature_job( new PoolMnPJob2 );
		mature_job->throw_ = throw_;
		mature_job->throw_bad_inputs_ = throw_bad_inputs_;
		mature_job->status_ = status_;

		core::pose::PoseOP pose = create_trpcage_ideal_poseop();
		mature_job->pose( pose );
		mature_job->mover( protocols::moves::NullMoverOP( new protocols::moves::NullMover ) );

		jobs_matured_.push_back( job->job_index() );

		return mature_job;
	}

	bool throw_;
	bool throw_bad_inputs_;
	JobStatus status_;

};


class PoolJQ3 : public PoolJQ1
{
public:
	PoolJQ3() :
		node_1_njobs_( 3 ),
		node_1_jobs_given_out_( false ),
		node_2_njobs_( 3 ),
		node_2_jobs_given_out_( false ),
		discarded_round_1_results_( false ),
		output_round_2_results_( false )
	{}

	JobDigraphOP initial_job_dag() override {
		JobDigraphOP job_dag( new JobDigraph( 2 ) );
		job_dag->add_edge( 1, 2 );
		return job_dag;
	}

	LarvalJobs determine_job_list( Size job_node_index, Size max_njobs ) override
	{
		LarvalJobs job_list;
		if ( job_node_index == 1 ) {
			if ( node_1_jobs_given_out_ ) return job_list;
			node_1_jobs_given_out_ = true;
			StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( node_1_njobs_, job_node_index ));
			job_list = expand_job_list( inner_job, max_njobs );
			//std::cout << "PoolJQ3 delivering round 1 jobs" << std::endl;
		} else if ( job_node_index == 2 ) {
			if ( node_2_jobs_given_out_ ) return job_list;
			node_2_jobs_given_out_ = true;
			StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( node_2_njobs_ ));
			inner_job->add_input_job_result_index( sp1(1) );
			job_list = expand_job_list( inner_job, max_njobs );

			//std::cout << "PoolJQ3 delivering round 2 jobs" << std::endl;
		}
		return job_list;
	}

	std::list< JobResultID > jobs_that_should_be_output() override
	{
		std::list< JobResultID > output_list;
		if ( node_2_jobs_given_out_ && ! output_round_2_results_ ) {
			bool all_jobs_completed = true;
			for ( Size ii = node_1_njobs_+1; ii <= node_1_njobs_ + node_2_njobs_; ++ii ) {
				if ( summaries_.count( sp1(ii) ) == 0 ) { all_jobs_completed = false; break;  }
			}
			if ( ! all_jobs_completed ) {
				return output_list;
			}
			for ( Size ii = node_1_njobs_+1; ii <= node_1_njobs_ + node_2_njobs_; ++ii ) {
				output_list.push_back( std::make_pair( ii, core::Size(1) ) );
			}
			output_round_2_results_ = true;
		}
		return output_list;
	}

	std::list< JobResultID > job_results_that_should_be_discarded() override
	{
		// After we have finished round 1, discard all the job results, except those
		// from job 1
		std::list< JobResultID > discard_list;

		if ( node_1_jobs_given_out_ && ! discarded_round_1_results_ ) {
			bool all_jobs_seen = true;
			for ( Size ii = 1; ii <= node_1_njobs_; ++ii ) {
				if ( summaries_.count( sp1(ii) ) == 0 ) { all_jobs_seen = false; break; }
			}
			if ( ! all_jobs_seen ) {
				return discard_list;
			}
			discarded_round_1_results_ = true;
			for ( Size ii = 2; ii <= node_1_njobs_; ++ii ) {
				discard_list.push_back( std::make_pair( ii, core::Size(1) ) );
			}
		}
		return discard_list;
	}



	Size node_1_njobs_;
	bool node_1_jobs_given_out_;

	Size node_2_njobs_;
	bool node_2_jobs_given_out_;

	bool discarded_round_1_results_;
	bool output_round_2_results_;

};

// The idea for this class is to test the various error messages delivered by the JD on node 0
// when a job result on an archive node can't be retrieved. In particular, if the job queen
// discards/outputs a job (removing it from the archive node) in between the time that the job
// distributor on node 0 sends the job out to a remote node to have it worked on and when that
// remote node tries to go and retrieve its results from the archive, then it's possible
// for the archive to legitimately not have the result.
class PoolJQ3b : public PoolJQ3
{
public:
	PoolJQ3b() : node_2_first_job_finished_( false ), discarded_job1_results_( false ), discard_job1_result_( true )
	{
		node_1_njobs_ = node_2_njobs_ = 2;
	}

	void note_job_completed( core::Size job_id, JobStatus status, Size nresults ) override {
		if ( job_id == node_1_njobs_ + 1 ) node_2_first_job_finished_ = true;
		PoolJQ1::note_job_completed( job_id, status, nresults );
	}

	std::list< JobResultID > jobs_that_should_be_output() override
	{
		if ( ! discard_job1_result_ && node_2_first_job_finished_ && ! discarded_job1_results_ ) {
			discarded_job1_results_ = true;
			std::list< JobResultID > job1_to_discard;
			job1_to_discard.push_back( sp( 1, 1 ) );
			return job1_to_discard;
		} else {
			return PoolJQ3::jobs_that_should_be_output();
		}
	}

	std::list< JobResultID > job_results_that_should_be_discarded() override
	{
		if ( discard_job1_result_ && node_2_first_job_finished_ && ! discarded_job1_results_ ) {
			discarded_job1_results_ = true;
			std::list< JobResultID > job1_to_discard;
			job1_to_discard.push_back( std::make_pair( 1, core::Size(1) ) );
			return job1_to_discard;
		} else {
			return PoolJQ3::job_results_that_should_be_discarded();
		}
	}

	bool node_2_first_job_finished_;
	bool discarded_job1_results_;
	bool discard_job1_result_; // either we ask the job distributor to discard job 1, or we ask the job distributor to output job 1
};

class PoolJQ4 : public PoolJQ1
{
public:
	using PoolJQ1::note_job_completed;
	using PoolJQ1::completed_job_summary;

	bool larval_job_needed_for_note_job_completed() const override { return true; }
	bool larval_job_needed_for_completed_job_summary() const override { return true; }

	void note_job_completed( LarvalJobCOP job, JobStatus status, Size nresults ) override {
		jobs_completed_through_larval_job_interface_.push_back( job->job_index() );
		PoolJQ1::note_job_completed( job->job_index(), status, nresults );
	}


	void completed_job_summary( LarvalJobCOP job, Size result_index, JobSummaryOP summary ) override {
		summaries_through_larval_job_interface_.push_back( job->job_index() );
		PoolJQ1::completed_job_summary( job->job_index(), result_index, summary );
	}


	utility::vector1< core::Size > jobs_completed_through_larval_job_interface_;
	utility::vector1< core::Size > summaries_through_larval_job_interface_;

};

class PoolJQ5 : public PoolJQ3
{
public:
	PoolJQ5() : node_3_njobs_( 1 ), node_3_jobs_given_out_( false ) { node_1_njobs_ = 1; node_2_njobs_ = 0; }

	JobDigraphOP initial_job_dag() override {
		JobDigraphOP job_dag( new JobDigraph( 1 ) );
		return job_dag;
	}

	void update_job_dag( JobDigraphUpdater & updater ) override {
		if ( updater.orig_num_nodes() == 1 && updater.job_dag()->get_job_node(1)->all_jobs_completed() ) {
			updater.add_node();
			updater.add_edge_to_new_node( 1, 2 );
		} else if ( updater.orig_num_nodes() == 2 && updater.job_dag()->get_job_node(2)->all_jobs_completed() ) {
			updater.add_node();
			updater.add_edge_to_new_node( 2, 3 );
		}
	}

	LarvalJobs determine_job_list( Size job_node_index, Size max_njobs ) override
	{
		if ( job_node_index < 3 ) return PoolJQ3::determine_job_list( job_node_index, max_njobs );
		LarvalJobs job_list;
		if ( ! node_3_jobs_given_out_ ) {
			node_3_jobs_given_out_ = true;
			StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( node_3_njobs_ ));
			job_list = expand_job_list( inner_job, max_njobs );
		}
		return job_list;
	}

	Size node_3_njobs_;
	bool node_3_jobs_given_out_;
};

// ok -- this class should be put into a larval job or
// a job result and when either the master node or the
// worker node goes to serialize it, serialization should
// fail. The JobDistributor should gracefully exit when this
// happens.
class Unserializable : public utility::pointer::ReferenceCount
{
public:
	Unserializable() : myint_( 5 ) {}
private:
	int myint_;
};

typedef utility::pointer::shared_ptr< Unserializable > UnserializableOP;

class Undeserializable : public utility::pointer::ReferenceCount
{
public:
	Undeserializable() : myint_( 5 ) {}
	int get_my_int() { return myint_; } // appease compiler
#ifdef SERIALIZATION

	template < class Archive >
	void
	save( Archive & arc ) const {
		arc( myint_ );
	}

	// this function causes trouble: it retrieves more data from the
	// archive than was put in
	template < class Archive >
	void
	load( Archive & ) {
		throw cereal::Exception( "Undeserializable could not be deserialized" );
	}

#endif

private:
	int myint_;
};

typedef utility::pointer::shared_ptr< Undeserializable > UndeserializableOP;


// This JQ puts an unserializable object into a LarvalJob.
class PoolJQ6 : public PoolJQ1
{
	// this class will imbed a non-serializable piece of data
	// in the larval job; the master node ought to catch an
	// exception when it tries to serialize the object and
	// send it to another node.
	LarvalJobs determine_job_list( Size job_node_index, Size max_njobs ) override
	{
		LarvalJobs jobs;
		if ( job_list_determined_ ) return jobs;

		StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( njobs_, job_node_index ) );
		inner_job->const_data_map().add( "testing", "testing", UnserializableOP( new Unserializable ));

		jobs = expand_job_list( inner_job, max_njobs );
		n_jobs_created_so_far_ += jobs.size();
		if ( n_jobs_created_so_far_ >= njobs_ ) {
			job_list_determined_ = true;
		}
		return jobs;
	}
};

// This JQ puts an undeserializable object into a LarvalJob.
class PoolJQ7 : public PoolJQ1
{
	// this class will imbed a non-serializable piece of data
	// in the larval job; the master node ought to catch an
	// exception when it tries to serialize the object and
	// send it to another node.
	LarvalJobs determine_job_list( Size, Size max_njobs ) override
	{
		LarvalJobs jobs;
		if ( job_list_determined_ ) return jobs;

		StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( njobs_ ) );
		inner_job->const_data_map().add( "testing", "testing", UndeserializableOP( new Undeserializable ));

		jobs = expand_job_list( inner_job, max_njobs );
		n_jobs_created_so_far_ += jobs.size();
		if ( n_jobs_created_so_far_ >= njobs_ ) {
			job_list_determined_ = true;
		}
		return jobs;
	}
};

// an unserializable job summary
class PoolJobSummary1 : public JobSummary
{
public:
	PoolJobSummary1() : JobSummary(), myint_( 5 ) {}
	int myint_;
};

// an undeserializable job summary
class PoolJobSummary2 : public JobSummary
{
public:
	PoolJobSummary2() : JobSummary(), myint_( 5 ) {}
	int get_my_int() { return myint_; } // appease compiler
#ifdef SERIALIZATION

	template < class Archive >
	void
	save( Archive & arc ) const {
		arc( myint_ );
	}

	// this function causes trouble: it retrieves more data from the
	// archive than was put in
	template < class Archive >
	void
	load( Archive & ) {
		throw cereal::Exception( "PoolJobSummary2 could not be deserialized" );
	}

#endif

	int myint_;
};

#ifdef SERIALIZATION

CEREAL_REGISTER_TYPE( Undeserializable )
CEREAL_REGISTER_TYPE( PoolJobSummary2 )

#endif

class PoolMnPJob3 : public MoverAndPoseJob
{
public:
	PoolMnPJob3() : MoverAndPoseJob() {}

	protocols::jd3::CompletedJobOutput
	run() override {
		CompletedJobOutput output = MoverAndPoseJob::run();
		PoseJobResultOP pose_result = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( output.job_results[1].second );
		pose_result->pose()->set_const_data( "testing", "testing", Unserializable() );
		return output;
	}
};

typedef utility::pointer::shared_ptr< PoolMnPJob3 > PoolMnPJob3OP;

class PoolMnPJob4 : public MoverAndPoseJob
{
public:
	PoolMnPJob4() : MoverAndPoseJob() {}

	protocols::jd3::CompletedJobOutput
	run() override {
		CompletedJobOutput output = MoverAndPoseJob::run();
		output.job_results[1].first = JobSummaryOP( new PoolJobSummary1 );
		return output;
	}
};

typedef utility::pointer::shared_ptr< PoolMnPJob4 > PoolMnPJob4OP;

class PoolMnPJob5 : public MoverAndPoseJob
{
public:
	PoolMnPJob5() : MoverAndPoseJob() {}

	CompletedJobOutput
	run() override {
		CompletedJobOutput output = MoverAndPoseJob::run();
		output.job_results[1].first = JobSummaryOP( new PoolJobSummary2 );
		return output;
	}
};

typedef utility::pointer::shared_ptr< PoolMnPJob5 > PoolMnPJob5OP;

class PoolJQ8 : public PoolJQ1
{

	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP job,
		utility::options::OptionCollectionCOP,
		utility::vector1< JobResultCOP > const &
	) override
	{
		using namespace protocols::jd3;
		PoolMnPJob3OP mature_job( new PoolMnPJob3 );

		core::pose::PoseOP pose = create_trpcage_ideal_poseop();
		mature_job->pose( pose );
		mature_job->mover( protocols::moves::NullMoverOP( new protocols::moves::NullMover ) );

		jobs_matured_.push_back( job->job_index() );

		return mature_job;
	}

};

class PoolJQ9 : public PoolJQ1
{

	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP job,
		utility::options::OptionCollectionCOP,
		utility::vector1< JobResultCOP > const &
	) override
	{
		using namespace protocols::jd3;
		PoolMnPJob4OP mature_job( new PoolMnPJob4 );

		core::pose::PoseOP pose = create_trpcage_ideal_poseop();
		mature_job->pose( pose );
		mature_job->mover( protocols::moves::NullMoverOP( new protocols::moves::NullMover ) );

		jobs_matured_.push_back( job->job_index() );

		return mature_job;
	}

};

class PoolJQ10 : public PoolJQ1
{

	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP job,
		utility::options::OptionCollectionCOP,
		utility::vector1< JobResultCOP > const &
	) override
	{
		using namespace protocols::jd3;
		PoolMnPJob5OP mature_job( new PoolMnPJob5 );

		core::pose::PoseOP pose = create_trpcage_ideal_poseop();
		mature_job->pose( pose );
		mature_job->mover( protocols::moves::NullMoverOP( new protocols::moves::NullMover ) );

		jobs_matured_.push_back( job->job_index() );

		return mature_job;
	}

};

class PoolJQ11 : public PoolJQ3
{
public:

	PoolJQ11() : distributed_deallocation_messages_( false ) {}

	std::list< deallocation::DeallocationMessageOP >
	deallocation_messages() override {
		// ok, let's send a message that two Poses should be deleted
		std::list< deallocation::DeallocationMessageOP > messages;
		if ( distributed_deallocation_messages_ ) return messages;
		distributed_deallocation_messages_ = true;

		typedef protocols::jd3::deallocation::InputPoseDeallocationMessage   PoseDeallocMsg;
		typedef protocols::jd3::deallocation::InputPoseDeallocationMessageOP PoseDeallocMsgOP;
		messages.push_back( PoseDeallocMsgOP( new PoseDeallocMsg( 1 ) ));
		messages.push_back( PoseDeallocMsgOP( new PoseDeallocMsg( 2 ) ));
		return messages;
	}

	void
	process_deallocation_message( deallocation::DeallocationMessageOP message ) override {
		received_messages_.push_back( message );
	}

	bool distributed_deallocation_messages_;
	utility::vector1< deallocation::DeallocationMessageOP > received_messages_;

};

typedef utility::pointer::shared_ptr< PoolJQ1 > PoolJQ1OP;
typedef utility::pointer::shared_ptr< PoolJQ2 > PoolJQ2OP;
typedef utility::pointer::shared_ptr< PoolJQ3 > PoolJQ3OP;
typedef utility::pointer::shared_ptr< PoolJQ3b > PoolJQ3bOP;
typedef utility::pointer::shared_ptr< PoolJQ4 > PoolJQ4OP;
typedef utility::pointer::shared_ptr< PoolJQ5 > PoolJQ5OP;
typedef utility::pointer::shared_ptr< PoolJQ6 > PoolJQ6OP;
typedef utility::pointer::shared_ptr< PoolJQ7 > PoolJQ7OP;
typedef utility::pointer::shared_ptr< PoolJQ8 > PoolJQ8OP;
typedef utility::pointer::shared_ptr< PoolJQ9 > PoolJQ9OP;
typedef utility::pointer::shared_ptr< PoolJQ10 > PoolJQ10OP;
typedef utility::pointer::shared_ptr< PoolJQ11 > PoolJQ11OP;



// --------------- Test Class --------------- //

class MPIWorkPoolJobDistributorTests : public CxxTest::TestSuite {

public:

	typedef core::Size   Size;
	typedef core::Real   Real;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
		basic::options::option[ basic::options::OptionKeys::jd3::compress_job_results ]( false );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

#ifdef SERIALIZATION

	template < class T >
	std::string
	serialized_T( T const & data ) {
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( data );
		}
		return oss.str();
	}


	template < class T >
	T
	deserialized_T( std::string const & str ) {
		std::istringstream iss( str );
		T return_object;
		cereal::BinaryInputArchive arc( iss );
		arc( return_object );
		return return_object;
	}

	LarvalJobOP
	deserialize_larval_job( std::string const & serialized_larval_job )
	{
		return deserialized_T< LarvalJobOP >( serialized_larval_job );
	}

	LarvalJobOP
	create_larval_job(
		core::Size njobs,
		core::Size nstruct_id,
		core::Size job_index
	)
	{
		StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( njobs ) );
		return LarvalJobOP( new LarvalJob( inner_job, nstruct_id, job_index ) );
	}

	// create a larval job that lists job 1 as a required input to this job
	LarvalJobOP
	create_larval_job_w_job1_dep(
		core::Size njobs,
		core::Size nstruct_id,
		core::Size job_index
	)
	{
		StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( njobs ) );
		inner_job->add_input_job_result_index( sp1(1) );
		LarvalJobOP larval_job( new LarvalJob( inner_job, nstruct_id, job_index ) );
		return larval_job;
	}


	std::string
	serialized_larval_job(
		core::Size njobs,
		core::Size nstruct_id,
		core::Size job_index
	)
	{
		LarvalJobOP larval_job = create_larval_job( njobs, nstruct_id, job_index );
		return serialized_T( larval_job );
	}

	std::string
	serialized_larval_job_and_job_result(
		core::Size njobs,
		core::Size nstruct_id,
		core::Size job_index,
		JobResultOP job_result
	)
	{
		return serialized_larval_job_and_job_result(
			create_larval_job( njobs, nstruct_id, job_index ), job_result );
	}

	std::string
	serialized_larval_job_and_job_result( LarvalJobOP larval_job, JobResultOP job_result )
	{
		return serialized_T( std::make_pair( larval_job, job_result ) );
	}

	std::string
	serialized_job_summaries( utility::vector1< JobSummaryOP > summaries )
	{
		return serialized_T( summaries );
	}

	utility::vector1< JobSummaryOP >
	deserialized_job_summaries( std::string const & str  )
	{
		return deserialized_T< utility::vector1< JobSummaryOP > >( str );
	}


	MPIWorkPoolJobDistributor::LarvalJobAndResult
	deserialized_larval_job_and_job_result( std::string const & str )
	{
		try {
			return deserialized_T< MPIWorkPoolJobDistributor::LarvalJobAndResult >( str );
		} catch ( ... ) {
			std::cerr << "Failed to deserialize larval job and job result!" << std::endl;
			TS_ASSERT( false );
			MPIWorkPoolJobDistributor::LarvalJobAndResult dummy_val;
			return dummy_val;
		}
	}

	std::list< deallocation::DeallocationMessageOP >
	deserialize_deallocation_msg_list( std::string const & serialized_msg_list )
	{
		return deserialized_T< std::list< deallocation::DeallocationMessageOP > >( serialized_msg_list );
	}


	std::string
	serialized_input_pose_deallocation_msg_list( std::list< core::Size > const & poses_to_deallocate )
	{
		typedef deallocation::DeallocationMessageOP MsgOP;
		typedef deallocation::InputPoseDeallocationMessage InputPoseMsg;
		std::list< MsgOP > msg_list;
		for ( auto pose_id : poses_to_deallocate ) {
			msg_list.push_back( MsgOP( new InputPoseMsg( pose_id ) ) );
		}
		return serialized_T( msg_list );
	}


#endif

	void test_MPIWorkPoolJobDistributor_master_node_vanilla_end_to_end_all_successes() {
		TS_ASSERT( true ); // appease non-serialization builds

#ifdef SERIALIZATION


		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		//CompletedJobOutput trpcage_job_output;
		//trpcage_job_output.first = trpcage_job_summary;
		//trpcage_job_output.second = trpcage_pose_result;

		//std::ostrinstream oss;
		//{
		// cereal::BinaryOutputArchive arc( oss );
		// arc( trpcage_job_output );
		//}
		//std::string serialized_trpcage_job_output = oss.str();

		SimulateMPI::initialize_simulation( 4 );

		// ok, all the nodes at once start requesting jobs from node 0
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );

		// let's start the job at node 1
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// and now the job at node 3
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// ok -- now let's pretend that node 2 has finished its job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );

		// ok -- before node 2 can finish its two part job-completion process,
		// node 1 will butt in and start its two part job-completion process,
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 2 ); // job_id
		send_size_to_node( 0, 1 ); // num jobs
		// ok -- node 0 is going to tell it to archive it's results on node 0

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 2 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 2, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num jobs

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 1, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 2 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 2 asks for a new job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 3, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 1 asks for a new job
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 3 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 3 asks for a new job
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// ok -- now we create a JQ on node 0, set it up to produce three jobs
		// and then create a job distributor for node 0 and tell it to go go go!

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ1OP jq( new PoolJQ1 );
		jq->njobs_ = 3;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 2 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv2 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized LarvalJob" );
		LarvalJobOP larval_job2 = deserialize_larval_job( ser_larv2 );
		TS_ASSERT_EQUALS( larval_job2->job_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 2 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 2 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 3 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv3 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized LarvalJob" );
		LarvalJobOP larval_job3 = deserialize_larval_job( ser_larv3 );
		TS_ASSERT_EQUALS( larval_job3->job_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 3 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 3 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq->status_.count( 2 ), 1 );
		TS_ASSERT_EQUALS( jq->status_.count( 3 ), 1 );
		TS_ASSERT_EQUALS( jq->status_[ 1 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq->status_[ 2 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq->status_[ 3 ], jd3_job_status_success );

		TS_ASSERT_EQUALS( jq->summaries_.count( sp1( 1 ) ), 1 );
		TS_ASSERT_EQUALS( jq->summaries_.count( sp1( 2 ) ), 1 );
		TS_ASSERT_EQUALS( jq->summaries_.count( sp1( 3 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq->summaries_[ sp1( 1 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq->summaries_[ sp1( 2 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq->summaries_[ sp1( 3 ) ])->energy(), 1234 );

		TS_ASSERT_EQUALS( jq->results_.count( sp1( 1 ) ), 1 );
		TS_ASSERT_EQUALS( jq->results_.count( sp1( 2 ) ), 1 );
		TS_ASSERT_EQUALS( jq->results_.count( sp1( 3 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq->results_[ sp1( 1 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq->results_[ sp1( 2 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq->results_[ sp1( 3 ) ])->pose()->total_residue(), 20 );

#endif
	}


	void test_workpool_jd3_worker_node_vanilla_end_to_end_all_successes() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank(0);

		// send 1st job
		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 2, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// tell node 1 where to send the results of job 1:
		send_integers_to_node( 1, utility::vector1< int >( 1, 0 ) );
		send_integer_to_node( 1, mpi_work_pool_jd_archival_completed );

		// Now send the second job to node 1
		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 2, 2, 2 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// tell node 1 where to send the results of job 2:
		send_integers_to_node( 1, utility::vector1< int >( 1, 0 ) );
		send_integer_to_node( 1, mpi_work_pool_jd_archival_completed );

		// Now tell node 1 to spin dow
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		try {
			SimulateMPI::set_mpi_rank( 1 );
			PoolJQ1OP jq1( new PoolJQ1 );
			MPIWorkPoolJobDistributor jd;
			jd.go( jq1 );
		} catch ( ... ) {
			std::cerr << "Exception thrown from worker node in a unit test that should not have thrown an exception!" << std::endl;
			TS_ASSERT( false );
			return;
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Then it replies that its job is complete
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 has finished its job", mpi_work_pool_jd_job_success );
		ts_assert_mpi_buffer_has_size( 1, "node 1 finished job 1", 1 );
		ts_assert_mpi_buffer_has_size( 1, "node 1 job 1 contained a single result", 1 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B2", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 wants to archive a job result", mpi_work_pool_jd_archive_job_result );
		ts_assert_mpi_buffer_has_size( 1, "node 1 wants to archive a result for job 1", 1 );
		ts_assert_mpi_buffer_has_size( 1, "node 1 wants to archive the first result for job 1", 1 );
		std::string serialized_job_result1 = ts_assert_mpi_buffer_has_string( 1, "node 1 sends the serialized trpcage pose result" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_and_result1 = deserialized_larval_job_and_job_result( serialized_job_result1 );
		TS_ASSERT( job_and_result1.first );
		TS_ASSERT( job_and_result1.second );
		if ( job_and_result1.first ) {
			TS_ASSERT_EQUALS( job_and_result1.first->job_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_max(), 2 );
		}
		PoseJobResultOP pose_result1 = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_and_result1.second );
		TS_ASSERT( pose_result1 );
		if ( pose_result1 ) {
			TS_ASSERT_EQUALS( pose_result1->pose()->total_residue(), 20 );
		}
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request C", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I finished archiving my completed job result -- I'm ready to send a JobSummary", mpi_work_pool_jd_job_success_and_archival_complete );
		ts_assert_mpi_buffer_has_size( 1, "node 1 says to node 0: this was for job index 1", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: job status is good", jd3_job_status_success );
		std::string serialized_summaries1 = ts_assert_mpi_buffer_has_string( 1, "node 1 says to node 0: here are the summaries" );
		utility::vector1< JobSummaryOP > summaries1 = deserialized_job_summaries( serialized_summaries1 );
		TS_ASSERT_EQUALS( summaries1.size(), 1 );
		if ( ! summaries1.size() == 1 ) return;
		TS_ASSERT( summaries1[1] );
		EnergyJobSummaryOP energy_summary1 = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( summaries1[1] );
		TS_ASSERT( energy_summary1 );
		if ( energy_summary1 ) {
			TS_ASSERT_EQUALS( energy_summary1->energy(), 1234 );
		}

		// Node 1 then should request a second job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request D", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Then it says the job is complete
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request E", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_job_success );
		ts_assert_mpi_buffer_has_size( 1, "node 1 finished job 2", 2 );
		ts_assert_mpi_buffer_has_size( 1, "node 1 job 2 had one result", 1 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request E2", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 needs to archive a result", mpi_work_pool_jd_archive_job_result );
		ts_assert_mpi_buffer_has_size( 1, "node 1 this is a result from job 2", 2 );
		ts_assert_mpi_buffer_has_size( 1, "node 1 this is the first result from job 2", 1 );
		std::string serialized_job_result2 = ts_assert_mpi_buffer_has_string( 1, "node 1 send the serialized trpcage pose result" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_and_result2 = deserialized_larval_job_and_job_result( serialized_job_result2 );
		TS_ASSERT( job_and_result2.first );
		TS_ASSERT( job_and_result2.second );
		if ( job_and_result2.first ) {
			TS_ASSERT_EQUALS( job_and_result2.first->job_index(), 2 );
			TS_ASSERT_EQUALS( job_and_result2.first->nstruct_index(), 2 );
			TS_ASSERT_EQUALS( job_and_result2.first->nstruct_max(), 2 );
		}
		PoseJobResultOP pose_result2 = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_and_result2.second );
		TS_ASSERT( pose_result2 );
		if ( pose_result2 ) {
			TS_ASSERT_EQUALS( pose_result2->pose()->total_residue(), 20 );
		}

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request F", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I finished archiving my completed job result -- I'm ready to send a JobSummary", mpi_work_pool_jd_job_success_and_archival_complete );
		ts_assert_mpi_buffer_has_size( 1, "node 1 says to node 0: this was for job index 2", 2 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that job 2 ended successfully", jd3_job_status_success );
		std::string serialized_summaries2 = ts_assert_mpi_buffer_has_string( 1, "node 1 says to node 0: here's the summary" );
		utility::vector1< JobSummaryOP > summaries2 = deserialized_job_summaries( serialized_summaries2 );
		TS_ASSERT_EQUALS( summaries2.size(), 1 );
		if ( summaries2.size() != 1 ) return;
		TS_ASSERT( summaries2[1] );
		EnergyJobSummaryOP energy_summary2 = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( summaries2[1] );
		TS_ASSERT( energy_summary2 );
		if ( energy_summary2 ) {
			TS_ASSERT_EQUALS( energy_summary2->energy(), 1234 );
		}

		// OK -- now node 1 should ask for a job for a final time
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// TO DO
		// make sure the job queen from node 1 has matured two larval jobs

#endif
	}

	void test_workpool_jd3_master_node_job_failed_do_not_retry() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_job_failed_do_not_retry ); // which is that the job I was running has failed
		send_size_to_node( 0, 1 ); // The job that failed was job #1

		// ask for another job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ1OP jq( new PoolJQ1 );
		jq->njobs_ = 1;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq->status_[ 1 ], jd3_job_status_failed_do_not_retry );
		TS_ASSERT_EQUALS( jq->summaries_.count( sp1(1) ), 0 );
		TS_ASSERT_EQUALS( jq->results_.count( sp1(1) ), 0 );

#endif

	}

	void test_workpool_jd3_master_node_job_failed_bad_input() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_job_failed_bad_input ); // which is that the job I was running has failed
		send_size_to_node( 0, 1 ); // The job that failed was job #1

		// ask for another job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ1OP jq( new PoolJQ1 );
		jq->njobs_ = 1;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq->status_[ 1 ], jd3_job_status_inputs_were_bad );
		TS_ASSERT_EQUALS( jq->summaries_.count( sp1(1) ), 0 );
		TS_ASSERT_EQUALS( jq->results_.count( sp1(1) ), 0 );

#endif
	}

	void test_workpool_jd3_master_node_job_failed_w_exception() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_job_failed_w_message ); // which is that the job I was running has failed
		send_size_to_node( 0, 1 ); // The job that failed was job #1
		send_string_to_node( 0, "This is a fake exception message" );

		// ask for another job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ1OP jq( new PoolJQ1 );
		jq->njobs_ = 1;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq->status_[ 1 ], jd3_job_status_failed_w_exception );
		TS_ASSERT_EQUALS( jq->summaries_.count( sp1(1) ), 0 );
		TS_ASSERT_EQUALS( jq->results_.count( sp1(1) ), 0 );

#endif
	}

	void test_workpool_jd3_master_node_job_failed_w_retry_limit_exceeded() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_job_failed_retry_limit_exceeded ); // which is that the job I was running has failed
		send_size_to_node( 0, 1 ); // The job that failed was job #1

		// ask for another job
		send_integer_to_node( 0, 1 ); // I have a request
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request ); // I want a new job

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ1OP jq( new PoolJQ1 );
		jq->njobs_ = 1;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq->status_[ 1 ], jd3_job_status_failed_max_retries );
		TS_ASSERT_EQUALS( jq->summaries_.count( sp1(1) ), 0 );
		TS_ASSERT_EQUALS( jq->results_.count( sp1(1) ), 0 );

#endif
	}

	void test_workpool_jd3_worker_node_job_failed_do_not_retry() {
		TS_ASSERT( true );


#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 0 );

		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// Now tell node 1 to spin dow
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		try {
			SimulateMPI::set_mpi_rank( 1 );
			PoolJQ2OP jq2( new PoolJQ2 );
			jq2->status_ = jd3_job_status_failed_do_not_retry;
			MPIWorkPoolJobDistributor jd;
			jd.go( jq2 );
		} catch ( ... ) {
			std::cerr << "Exception thrown from worker node in a unit test that should not have thrown an exception!" << std::endl;
			TS_ASSERT( false );
			return;
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Now node 1 should have said that the job failed
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that the job failed with a status 'do not retry'", mpi_work_pool_jd_job_failed_do_not_retry );
		ts_assert_mpi_buffer_has_size( 1, "node 1 says that the job which failed had index 1", 1 );

		// Node 1 requests another job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request C", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests another job of node 0", mpi_work_pool_jd_new_job_request );
#endif
	}

	void test_workpool_jd3_worker_node_job_failed_bad_input() {
		TS_ASSERT( true );


#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 0 );

		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// Now tell node 1 to spin dow
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		try {
			SimulateMPI::set_mpi_rank( 1 );
			PoolJQ2OP jq2( new PoolJQ2 );
			jq2->status_ = jd3_job_status_inputs_were_bad;
			MPIWorkPoolJobDistributor jd;
			jd.go( jq2 );
		} catch ( ... ) {
			std::cerr << "Exception thrown from worker node in a unit test that should not have thrown an exception!" << std::endl;
			TS_ASSERT( false );
			return;
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Now node 1 should have said that the job failed
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that the job failed with a status 'do not retry'", mpi_work_pool_jd_job_failed_bad_input );
		ts_assert_mpi_buffer_has_size( 1, "node 1 says that the job which failed had index 1", 1 );

		// Node 1 requests another job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request C", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests another job of node 0", mpi_work_pool_jd_new_job_request );
#endif
	}

	void test_workpool_jd3_worker_node_job_failed_w_exception() {
		TS_ASSERT( true );


#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 0 );

		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// Now tell node 1 to spin dow
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		try {
			SimulateMPI::set_mpi_rank( 1 );
			PoolJQ2OP jq2( new PoolJQ2 );
			jq2->throw_ = true;
			MPIWorkPoolJobDistributor jd;
			jd.go( jq2 );
		} catch ( ... ) {
			std::cerr << "Exception thrown from worker node in a unit test that should not have thrown an exception!" << std::endl;
			TS_ASSERT( false );
			return;
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Now node 1 should have said that the job failed
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that the job failed with exception message'", mpi_work_pool_jd_job_failed_w_message );
		ts_assert_mpi_buffer_has_size( 1, "node 1 says that the job which failed had index 1", 1 );
		ts_assert_mpi_buffer_has_string( 1, "node 1 hands the exception message to node 0", "PoolMnPJob2 exception" );

		// Node 1 requests another job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request C", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests another job of node 0", mpi_work_pool_jd_new_job_request );
#endif
	}

	void test_workpool_jd3_worker_node_job_failed_w_bad_inputs_exception() {
		TS_ASSERT( true );


#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 0 );

		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// Now tell node 1 to spin dow
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		try {
			SimulateMPI::set_mpi_rank( 1 );
			PoolJQ2OP jq2( new PoolJQ2 );
			jq2->throw_bad_inputs_ = true;
			MPIWorkPoolJobDistributor jd;
			jd.go( jq2 );
		} catch ( ... ) {
			std::cerr << "Exception thrown from worker node in a unit test that should not have thrown an exception!" << std::endl;
			TS_ASSERT( false );
			return;
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Now node 1 should have said that the job failed
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that the job failed with exception message'", mpi_work_pool_jd_job_failed_w_message );
		ts_assert_mpi_buffer_has_size( 1, "node 1 says that the job which failed had index 1", 1 );
		ts_assert_mpi_buffer_has_string( 1, "node 1 hands the exception message to node 0", "PoolMnPJob2 bad input exception" );

		// Node 1 requests another job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request C", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests another job of node 0", mpi_work_pool_jd_new_job_request );
#endif
	}

	void test_workpool_jd3_worker_node_job_failed_w_retry_limit_exceeded() {
		TS_ASSERT( true );


#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 0 );

		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		LarvalJobOP multi_pass_larval_job = create_larval_job( 1, 1, 1 );
		multi_pass_larval_job->retry_limit( 5 ); // try this job multiple times
		send_string_to_node( 1, serialized_T( multi_pass_larval_job ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// Now tell node 1 to spin dow
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		try {
			SimulateMPI::set_mpi_rank( 1 );
			PoolJQ2OP jq2( new PoolJQ2 );
			jq2->status_ = jd3_job_status_failed_retry;
			MPIWorkPoolJobDistributor jd;
			jd.go( jq2 );
		} catch ( ... ) {
			std::cerr << "Exception thrown from worker node in a unit test that should not have thrown an exception!" << std::endl;
			TS_ASSERT( false );
			return;
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Now node 1 should have said that the job failed
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that the job failed after hitting its retry limit'", mpi_work_pool_jd_job_failed_retry_limit_exceeded );
		ts_assert_mpi_buffer_has_size( 1, "node 1 says that the job which failed had index 1", 1 );

		// Node 1 requests another job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request C", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests another job of node 0", mpi_work_pool_jd_new_job_request );
#endif
	}

	void test_jd3_workpool_manager_two_rounds_no_archive() {


#ifdef SERIALIZATION

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 4 );

		// ok, all the nodes at once start requesting jobs from node 0
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );

		// let's start the job at node 1
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// and now the job at node 3
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// ok -- now let's pretend that node 2 has finished its job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );

		// ok -- before node 2 can finish its two part job-completion process,
		// node 1 will butt in and start its two part job-completion process,
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 2 ); // job_id
		send_size_to_node( 0, 1 ); // num jobs
		// ok -- node 0 is going to tell it to archive it's results on node 0

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 2 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 2, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num jobs

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 1, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 2 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // num results


		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 3, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 1 asks for a new job
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 2 );
		// Now node 2 asks for a new job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );


		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 3 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 3 asks for a new job
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 1 );

		// now retrieve the job result for job 1
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 0, 1 ); // retrieve the result from job #1
		send_size_to_node( 0, 1 ); // result #1 from job #1

		SimulateMPI::set_mpi_rank( 2 );

		// now retrieve the job result for job 1
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 0, 1 ); // retrieve the result from job #1
		send_size_to_node( 0, 1 ); // result #1 from job #1

		SimulateMPI::set_mpi_rank( 3 );

		// now retrieve the job result for job 1
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 0, 1 ); // retrieve the result from job #1
		send_size_to_node( 0, 1 ); // result #1 from job #1

		// Now the second round of job results start trickling in

		SimulateMPI::set_mpi_rank( 1 );

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 4 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 4 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 4, trpcage_pose_result ));

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 4 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 6 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 6 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 6, trpcage_pose_result ));

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 6 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 5 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 5 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 5, trpcage_pose_result ));

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 5 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		//OK! Now create a PoolJQ3 and let 'er rip
		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ3OP jq3( new PoolJQ3 );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3 );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv4 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job4 = deserialize_larval_job( ser_larv4 );
		TS_ASSERT_EQUALS( larval_job4->job_index(), 4 );
		TS_ASSERT_EQUALS( larval_job4->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job4->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 a vector1 of job result indices with one element whose value is 0", utility::vector1< Size >( 1, 0 ) );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1a = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1a = deserialized_larval_job_and_job_result( ser_job_res1a );
		TS_ASSERT( job_res1a.first );
		TS_ASSERT( job_res1a.second );
		TS_ASSERT_EQUALS( job_res1a.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1a.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1a.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1a.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 2 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv2 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized LarvalJob" );
		LarvalJobOP larval_job2 = deserialize_larval_job( ser_larv2 );
		TS_ASSERT_EQUALS( larval_job2->job_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 2 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 2 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv5 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized LarvalJob" );
		LarvalJobOP larval_job5 = deserialize_larval_job( ser_larv5 );
		TS_ASSERT_EQUALS( larval_job5->job_index(), 5 );
		TS_ASSERT_EQUALS( larval_job5->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( larval_job5->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 2 a vector1 of job result indices with one element whose value is 0", utility::vector1< Size >( 1, 0 ) );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1b = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1b = deserialized_larval_job_and_job_result( ser_job_res1b );
		TS_ASSERT( job_res1b.first );
		TS_ASSERT( job_res1b.second );
		TS_ASSERT_EQUALS( job_res1b.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1b.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1b.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1b.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 2 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 3 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv3 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized LarvalJob" );
		LarvalJobOP larval_job3 = deserialize_larval_job( ser_larv3 );
		TS_ASSERT_EQUALS( larval_job3->job_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 3 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 3 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv6 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized LarvalJob" );
		LarvalJobOP larval_job6 = deserialize_larval_job( ser_larv6 );
		TS_ASSERT_EQUALS( larval_job6->job_index(), 6 );
		TS_ASSERT_EQUALS( larval_job6->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( larval_job6->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 3 a vector1 of job result indices with one element whose value is 0", utility::vector1< Size >( 1, 0 ) );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1c = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1c = deserialized_larval_job_and_job_result( ser_job_res1b );
		TS_ASSERT( job_res1c.first );
		TS_ASSERT( job_res1c.second );
		TS_ASSERT_EQUALS( job_res1c.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1c.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1c.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1c.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 3 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq3->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 2 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 3 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 4 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 5 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 6 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_[ 1 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 2 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 3 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 4 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 5 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 6 ], jd3_job_status_success );

		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 1 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 2 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 3 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 4 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 5 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 6 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 1 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 2 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 3 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 4 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 5 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 6 ) ])->energy(), 1234 );

		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 1 ) ), 0 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 2 ) ), 0 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 3 ) ), 0 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 4 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 5 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 6 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq3->results_[ sp1( 4 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq3->results_[ sp1( 5 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq3->results_[ sp1( 6 ) ])->pose()->total_residue(), 20 );



#endif
	}

	void test_jd3_workpool_manager_two_rounds_one_archive_no_archiving_on_node_0() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );


		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 5 );

		// ok, all the nodes at once start requesting jobs from node 0
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );
		SimulateMPI::set_mpi_rank( 4 );
		send_integer_to_node( 0, 4 );

		// let's start the job at node 2
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// and now the job at node 4
		SimulateMPI::set_mpi_rank( 4 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// ok -- now let's pretend that node 3 has finished its job
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );

		// ok -- before node 3 can finish its two part job-completion process,
		// node 2 will butt in and start its two part job-completion process,
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 2 ); // job_id
		send_size_to_node( 0, 1 ); // num results
		// sent to archive instead -- send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 2, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// sent to archive instead -- send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 1, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 2 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 4 );
		send_integer_to_node( 0, 4 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// sent to archive instead -- send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 3, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // job id
		send_integer_to_node( 0, jd3_job_status_success ); // job id
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 2 asks for a new job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		// Now node 3 asks for a new job
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );


		SimulateMPI::set_mpi_rank( 4 );
		send_integer_to_node( 0, 4 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 3 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 4 asks for a new job
		send_integer_to_node( 0, 4 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 2 );

		// now retrieve the job result for job 1
		// send_integer_to_node( 0, 2 );
		// send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		// send_size_to_node( 0, 1 ); // retrieve the result from job #1

		SimulateMPI::set_mpi_rank( 3 );

		// now retrieve the job result for job 1
		// send_integer_to_node( 0, 3 );
		// send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		// send_size_to_node( 0, 1 ); // retrieve the result from job #1

		SimulateMPI::set_mpi_rank( 4 );

		// now retrieve the job result for job 1
		// send_integer_to_node( 0, 4 );
		// send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		// send_size_to_node( 0, 1 ); // retrieve the result from job #1

		// Now the second round of job results start trickling in

		SimulateMPI::set_mpi_rank( 2 );

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 4 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// sent to archive instead -- send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 4, trpcage_pose_result ));

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 4 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 4 );
		send_integer_to_node( 0, 4 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 6 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// sent to archive instead -- send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 6, trpcage_pose_result ));

		send_integer_to_node( 0, 4 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 6 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 4 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 5 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// sent to archive instead -- send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 5, trpcage_pose_result ));

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 5 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 1 );
		// ok, now the archive is going to receive a series of job result requests from
		// node 0 and respond by sending back the serialized larval-jobs and job results
		send_integer_to_node( 0, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 4, trpcage_pose_result ));

		send_integer_to_node( 0, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 5, trpcage_pose_result ));

		send_integer_to_node( 0, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 6, trpcage_pose_result ));

		//OK! Now create a PoolJQ3 and let 'er rip
		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ3OP jq3( new PoolJQ3 );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3 );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 2 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 2 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 2 to archive its results on node 1", utility::vector1< int >( 1, 1 ) );
		// ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv4 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized LarvalJob" );
		LarvalJobOP larval_job4 = deserialize_larval_job( ser_larv4 );
		TS_ASSERT_EQUALS( larval_job4->job_index(), 4 );
		TS_ASSERT_EQUALS( larval_job4->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job4->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 2 a vector1 of job result indices with one element whose value is 1", utility::vector1< Size >( 1, 1 ) );

		//ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		// std::string ser_job_res1a = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		// MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1a = deserialized_larval_job_and_job_result( ser_job_res1a );
		// TS_ASSERT( job_res1a.first );
		// TS_ASSERT( job_res1a.second );
		// TS_ASSERT_EQUALS( job_res1a.first->job_index(), 1 );
		// TS_ASSERT_EQUALS( job_res1a.first->nstruct_index(), 1 );
		// TS_ASSERT_EQUALS( job_res1a.first->nstruct_max(), 3 );
		// TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1a.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 2 to archive its results on node 1", utility::vector1< int >( 1, 1 ) );
		// ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 3 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv2 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized LarvalJob" );
		LarvalJobOP larval_job2 = deserialize_larval_job( ser_larv2 );
		TS_ASSERT_EQUALS( larval_job2->job_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 3 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 3 to archive its results on node 1", utility::vector1< int >( 1, 1 ) );
		//ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv5 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized LarvalJob" );
		LarvalJobOP larval_job5 = deserialize_larval_job( ser_larv5 );
		TS_ASSERT_EQUALS( larval_job5->job_index(), 5 );
		TS_ASSERT_EQUALS( larval_job5->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( larval_job5->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 3 a vector1 of job result indices with one element whose value is 1", utility::vector1< Size >( 1, 1 ) );

		// ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		// std::string ser_job_res1b = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		// MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1b = deserialized_larval_job_and_job_result( ser_job_res1b );
		// TS_ASSERT( job_res1b.first );
		// TS_ASSERT( job_res1b.second );
		// TS_ASSERT_EQUALS( job_res1b.first->job_index(), 1 );
		// TS_ASSERT_EQUALS( job_res1b.first->nstruct_index(), 1 );
		// TS_ASSERT_EQUALS( job_res1b.first->nstruct_max(), 3 );
		// TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1b.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 3 to archive its results on node 1", utility::vector1< int >( 1, 1 ) );
		// ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 4 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 4 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv3 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 4 the serialized LarvalJob" );
		LarvalJobOP larval_job3 = deserialize_larval_job( ser_larv3 );
		TS_ASSERT_EQUALS( larval_job3->job_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 4 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 4 to archive its results on node 1", utility::vector1< int >( 1, 1 ) );
		// ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 4 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 4 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv6 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 4 the serialized LarvalJob" );
		LarvalJobOP larval_job6 = deserialize_larval_job( ser_larv6 );
		TS_ASSERT_EQUALS( larval_job6->job_index(), 6 );
		TS_ASSERT_EQUALS( larval_job6->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( larval_job6->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 4 a vector1 of job result indices with one element whose value is 1", utility::vector1< Size >( 1, 1 ) );

		// ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 4 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		// std::string ser_job_res1c = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		// MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1c = deserialized_larval_job_and_job_result( ser_job_res1b );
		// TS_ASSERT( job_res1c.first );
		// TS_ASSERT( job_res1c.second );
		// TS_ASSERT_EQUALS( job_res1c.first->job_index(), 1 );
		// TS_ASSERT_EQUALS( job_res1c.first->nstruct_index(), 1 );
		// TS_ASSERT_EQUALS( job_res1c.first->nstruct_max(), 3 );
		// TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1c.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 4 to archive its results on node 1", utility::vector1< int >( 1, 1 ) );
		// ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 4 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 4 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 1 ); // the archive node
		// first messages from node 0 are to discard the results from jobs 2 and 3
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it wants its attention A", 0 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it needs to discard a job A", mpi_work_pool_jd_discard_job_result );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 2", 2 );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 2 result 1", 1 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it wants its attention B", 0 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it needs to discard a job B", mpi_work_pool_jd_discard_job_result );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 3", 3 );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 3 result 1", 1 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it wants its attention C", 0 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it needs to retrieve a job A", mpi_work_pool_jd_retrieve_and_discard_job_result );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 4", 4 );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 4 result 1", 1 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it wants its attention D", 0 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it needs to retrieve a job B", mpi_work_pool_jd_retrieve_and_discard_job_result );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 5", 5 );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 5 result 1", 1 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it wants its attention E", 0 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it needs to retrieve a job C", mpi_work_pool_jd_retrieve_and_discard_job_result );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 6", 6 );
		ts_assert_mpi_buffer_has_size( 0, "node 0 tells node 1 to discard job 6 result 1", 1 );

		TS_ASSERT_EQUALS( jq3->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 2 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 3 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 4 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 5 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_.count( 6 ), 1 );
		TS_ASSERT_EQUALS( jq3->status_[ 1 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 2 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 3 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 4 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 5 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq3->status_[ 6 ], jd3_job_status_success );

		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 1 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 2 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 3 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 4 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 5 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->summaries_.count( sp1( 6 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 1 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 2 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 3 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 4 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 5 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq3->summaries_[ sp1( 6 ) ])->energy(), 1234 );

		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 1 ) ), 0 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 2 ) ), 0 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 3 ) ), 0 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 4 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 5 ) ), 1 );
		TS_ASSERT_EQUALS( jq3->results_.count( sp1( 6 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq3->results_[ sp1( 4 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq3->results_[ sp1( 5 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq3->results_[ sp1( 6 ) ])->pose()->total_residue(), 20 );

#endif
	}

	void test_jd3_workpool_archive_two_rounds_one_archive_no_archiving_on_node_0() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 5 );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 1, 3 ); // node 3 says "I have a message for you, archive"
		send_integer_to_node( 1, mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		send_size_to_node( 1, 2 ); // job_id
		send_size_to_node( 1, 1 ); // result index
		send_string_to_node( 1, serialized_larval_job_and_job_result( 3, 2, 2, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 1, 2 ); // node 2 says "I have a message for you, archive"
		send_integer_to_node( 1, mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		send_size_to_node( 1, 1 ); // job_id
		send_size_to_node( 1, 1 ); // result index
		send_string_to_node( 1, serialized_larval_job_and_job_result( 3, 1, 1, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 4 );
		send_integer_to_node( 1, 4 ); // node 4 says "I have a message for you, archive"
		send_integer_to_node( 1, mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		send_size_to_node( 1, 3 ); // job_id
		send_size_to_node( 1, 1 ); // result index
		send_string_to_node( 1, serialized_larval_job_and_job_result( 3, 3, 3, trpcage_pose_result ));


		SimulateMPI::set_mpi_rank( 0 );
		// master node says "discard job results 2 and 3
		send_integer_to_node( 1, 0 ); // node 0 says to archive node: "hey buddy, we need to talk"
		send_integer_to_node( 1, mpi_work_pool_jd_discard_job_result ); // "I need you to discard a job result"
		send_size_to_node( 1, 2 ); // discard job 2
		send_size_to_node( 1, 1 ); // discard job 2 result index 1

		send_integer_to_node( 1, 0 ); // node 0 says to archive node: "hey buddy, we need to talk"
		send_integer_to_node( 1, mpi_work_pool_jd_discard_job_result ); // "I need you to discard a job result"
		send_size_to_node( 1, 3 ); // discard job 3
		send_size_to_node( 1, 1 ); // discard job 3 result index 1

		SimulateMPI::set_mpi_rank( 2 );

		// now retrieve the job result for job 1
		send_integer_to_node( 1, 2 );
		send_integer_to_node( 1, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 1, 1 ); // retrieve the result from job #1
		send_size_to_node( 1, 1 ); // retrieve the result from job #1 result index #1

		SimulateMPI::set_mpi_rank( 3 );

		// now retrieve the job result for job 1
		send_integer_to_node( 1, 3 );
		send_integer_to_node( 1, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 1, 1 ); // retrieve the result from job #1
		send_size_to_node( 1, 1 ); // retrieve the result from job #1 result index #1

		SimulateMPI::set_mpi_rank( 4 );

		// now retrieve the job result for job 1
		send_integer_to_node( 1, 4 );
		send_integer_to_node( 1, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 1, 1 ); // retrieve the result from job #1
		send_size_to_node( 1, 1 ); // retrieve the result from job #1 result index #1

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 1, 2 ); // node 2 says "I have a message for you, archive"
		send_integer_to_node( 1, mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		send_size_to_node( 1, 4 ); // job_id
		send_size_to_node( 1, 1 ); // result index #1
		send_string_to_node( 1, serialized_larval_job_and_job_result( 3, 1, 4, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 1, 3 ); // node 3 says "I have a message for you, archive"
		send_integer_to_node( 1, mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		send_size_to_node( 1, 5 ); // job_id
		send_size_to_node( 1, 1 ); // result index #1
		send_string_to_node( 1, serialized_larval_job_and_job_result( 3, 2, 5, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 4 );
		send_integer_to_node( 1, 4 ); // node 4 says "I have a message for you, archive"
		send_integer_to_node( 1, mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		send_size_to_node( 1, 6 ); // job_id
		send_size_to_node( 1, 1 ); // result index #1
		send_string_to_node( 1, serialized_larval_job_and_job_result( 3, 3, 6, trpcage_pose_result ));

		// Finally, node 0 asks the archive for job result

		SimulateMPI::set_mpi_rank( 0 );

		send_integer_to_node( 1, 0 ); // node 0 says "hey, I need to talk to you"
		send_integer_to_node( 1, mpi_work_pool_jd_retrieve_and_discard_job_result ); // I need a job result
		send_size_to_node( 1, 4 ); // job id
		send_size_to_node( 1, 1 ); // result index

		send_integer_to_node( 1, 0 ); // node 0 says "hey, I need to talk to you"
		send_integer_to_node( 1, mpi_work_pool_jd_retrieve_and_discard_job_result ); // I need a job result
		send_size_to_node( 1, 5 ); // job id
		send_size_to_node( 1, 1 ); // result index

		send_integer_to_node( 1, 0 ); // node 0 says "hey, I need to talk to you"
		send_integer_to_node( 1, mpi_work_pool_jd_retrieve_and_discard_job_result ); // I need a job result
		send_size_to_node( 1, 6 ); // job id
		send_size_to_node( 1, 1 ); // result index

		send_integer_to_node( 1, 0 ); // node 0 says "hey, I need to talk to you"
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down ); // time to shut down

		//OK! Now create a PoolJQ3 and let 'er rip
		SimulateMPI::set_mpi_rank( 1 );
		PoolJQ3OP jq3( new PoolJQ3 );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3 );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 1

		SimulateMPI::set_mpi_rank( 2 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 2 that archival completed (for job 2)", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 2 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1a = ts_assert_mpi_buffer_has_string( 1, "node 1 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1a = deserialized_larval_job_and_job_result( ser_job_res1a );
		TS_ASSERT( job_res1a.first );
		TS_ASSERT( job_res1a.second );
		TS_ASSERT_EQUALS( job_res1a.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1a.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1a.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1a.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 2 that archival completed (for job 4)", mpi_work_pool_jd_archival_completed );

		SimulateMPI::set_mpi_rank( 3 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 3 that archival completed (for job 1)", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 3 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1b = ts_assert_mpi_buffer_has_string( 1, "node 1 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1b = deserialized_larval_job_and_job_result( ser_job_res1b );
		TS_ASSERT( job_res1b.first );
		TS_ASSERT( job_res1b.second );
		TS_ASSERT_EQUALS( job_res1b.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1b.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1b.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1b.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 3 that archival completed (for job 5)", mpi_work_pool_jd_archival_completed );

		SimulateMPI::set_mpi_rank( 4 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 4 that archival completed (for job 3)", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 4 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1c = ts_assert_mpi_buffer_has_string( 1, "node 1 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1c = deserialized_larval_job_and_job_result( ser_job_res1b );
		TS_ASSERT( job_res1c.first );
		TS_ASSERT( job_res1c.second );
		TS_ASSERT_EQUALS( job_res1c.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1c.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1c.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1c.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 4 that archival completed", mpi_work_pool_jd_archival_completed );

		SimulateMPI::set_mpi_rank( 0 ); // the master node

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that it has job result 4", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res4 = ts_assert_mpi_buffer_has_string( 1, "node 1 sends the serialized larval-job and job-result from job 4" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res4 = deserialized_larval_job_and_job_result( ser_job_res4 );
		TS_ASSERT( job_res4.first );
		TS_ASSERT( job_res4.second );
		TS_ASSERT_EQUALS( job_res4.first->job_index(), 4 );
		TS_ASSERT_EQUALS( job_res4.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res4.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1a.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that it has job result 5", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res5 = ts_assert_mpi_buffer_has_string( 1, "node 1 sends the serialized larval-job and job-result from job 5" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res5 = deserialized_larval_job_and_job_result( ser_job_res5 );
		TS_ASSERT( job_res5.first );
		TS_ASSERT( job_res5.second );
		TS_ASSERT_EQUALS( job_res5.first->job_index(), 5 );
		TS_ASSERT_EQUALS( job_res5.first->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( job_res5.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1a.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that it has job result 6", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res6 = ts_assert_mpi_buffer_has_string( 1, "node 1 sends the serialized larval-job and job-result from job 6" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res6 = deserialized_larval_job_and_job_result( ser_job_res6 );
		TS_ASSERT( job_res6.first );
		TS_ASSERT( job_res6.second );
		TS_ASSERT_EQUALS( job_res6.first->job_index(), 6 );
		TS_ASSERT_EQUALS( job_res6.first->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( job_res6.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1a.second )->pose()->total_residue(), 20 );


#endif
	}

	void test_jd3_workpool_worker_two_rounds_one_archive_no_archiving_on_node_0() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 3 );

		SimulateMPI::set_mpi_rank( 0 );

		// send 1st job
		send_integer_to_node( 2, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 2, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 2, utility::vector1< Size >() );

		// tell node 2 where to send the results of job 1:
		send_integers_to_node( 2, utility::vector1< int > (1,1) );

		SimulateMPI::set_mpi_rank( 1 );
		// ok, node 2 will say that it needs to send us a job, and it'll send it to node 1 (the archive)
		// and then the archive needs to reply that the archival was successful
		send_integer_to_node( 2, mpi_work_pool_jd_archival_completed );

		SimulateMPI::set_mpi_rank( 0 );
		// send 2nd job, which lists job 1 as a required input
		send_integer_to_node( 2, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 2, serialized_T( create_larval_job_w_job1_dep( 1, 1, 2 ) ));
		send_sizes_to_node( 2, utility::vector1< Size >( 1, 1 ) ); // say that job 1's results live on node 1

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 2, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( 2, serialized_larval_job_and_job_result( 1, 1, 1, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 0 );
		// tell node 2 where to send the results of job 2:
		send_integers_to_node( 2, utility::vector1< int > (1,1) );

		SimulateMPI::set_mpi_rank( 1 );
		// ok, node 2 will say that it needs to send us a job, and it'll send it to node 1 (the archive)
		// and then the archive needs to reply that the archival was successful
		send_integer_to_node( 2, mpi_work_pool_jd_archival_completed );

		SimulateMPI::set_mpi_rank( 0 );
		send_integer_to_node( 2, mpi_work_pool_jd_spin_down ); // no more jobs to hand out

		//OK! Now create a PoolJQ3 and let 'er rip
		SimulateMPI::set_mpi_rank( 2 );
		PoolJQ3OP jq3( new PoolJQ3 );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3 );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 2

		SimulateMPI::set_mpi_rank( 0 );
		// Node 2 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request A", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Then it replies that its job is complete
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request B", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 has finished its job", mpi_work_pool_jd_job_success );
		ts_assert_mpi_buffer_has_size( 2, "node 2 finished job #1", 1 );
		ts_assert_mpi_buffer_has_size( 2, "job #1 generated 1 result", 1 );

		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request C", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I finished archiving my completed job result -- I'm ready to send a JobSummary", mpi_work_pool_jd_job_success_and_archival_complete );
		ts_assert_mpi_buffer_has_size( 2, "node 2 says to node 0: this was for job index 1", 1 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 reports the job status to node 0: success", jd3_job_status_success );
		std::string serialized_summaries1 = ts_assert_mpi_buffer_has_string( 2, "node 2 says to node 0: here's the summary" );
		utility::vector1< JobSummaryOP > summaries1 = deserialized_job_summaries( serialized_summaries1 );
		TS_ASSERT_EQUALS( summaries1.size(), 1 );
		if ( summaries1.size() != 1 ) return;
		TS_ASSERT( summaries1[1] );
		EnergyJobSummaryOP energy_summary1 = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( summaries1[1] );
		TS_ASSERT( energy_summary1 );
		if ( energy_summary1 ) {
			TS_ASSERT_EQUALS( energy_summary1->energy(), 1234 );
		}

		// Now Node 2 requests another job
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request D", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Then it replies that its job is complete
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request E", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 has finished its job", mpi_work_pool_jd_job_success );
		ts_assert_mpi_buffer_has_size( 2, "node 2 finished job #2", 2 );
		ts_assert_mpi_buffer_has_size( 2, "job #2 produced one result", 1 );

		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request F", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I finished archiving my completed job result -- I'm ready to send a JobSummary", mpi_work_pool_jd_job_success_and_archival_complete );
		ts_assert_mpi_buffer_has_size( 2, "node 2 says to node 0: this was for job index 2", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "job 2 was a success", jd3_job_status_success );
		std::string serialized_summaries2 = ts_assert_mpi_buffer_has_string( 2, "node 2 says to node 0: here are the summaries" );
		utility::vector1< JobSummaryOP > summaries2 = deserialized_job_summaries( serialized_summaries2 );
		TS_ASSERT_EQUALS( summaries2.size(), 1 );
		if ( summaries2.size() != 1 ) return;
		TS_ASSERT( summaries2[ 1 ] );
		EnergyJobSummaryOP energy_summary2 = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( summaries2[1] );
		TS_ASSERT( energy_summary2 );
		if ( energy_summary2 ) {
			TS_ASSERT_EQUALS( energy_summary2->energy(), 1234 );
		}

		// Now Node 2 requests another job, and it will turn out that there are no more, so node 0 will send a spin down signal.
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request G", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 1 );

		// Let's see what node 2 sent to the archive node.
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"I have a message for you, archive\"", 2 ); // node 2 says "I have a message for you, archive"
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"please archive this job result\"", mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		ts_assert_mpi_buffer_has_size( 2, "node 2 says this job result is for job #1", 1 ); // job_id
		ts_assert_mpi_buffer_has_size( 2, "node 2 says this job result is for job #1 result #1", 1 ); // result index
		std::string serialized_job_result1 = ts_assert_mpi_buffer_has_string( 2, "node 2 sends the serialized trpcage pose result for job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_and_result1 = deserialized_larval_job_and_job_result( serialized_job_result1 );
		TS_ASSERT( job_and_result1.first );
		TS_ASSERT( job_and_result1.second );
		if ( job_and_result1.first ) {
			TS_ASSERT_EQUALS( job_and_result1.first->job_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_max(), 1 );
		}
		PoseJobResultOP pose_result1 = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_and_result1.second );
		TS_ASSERT( pose_result1 );
		if ( pose_result1 ) {
			TS_ASSERT_EQUALS( pose_result1->pose()->total_residue(), 20 );
		}

		// node 2 then asks for the job results from job 1
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to the archive that it has a request", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says that it needs a job result",  mpi_work_pool_jd_retrieve_job_result );
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that it wants the job with index 1", 1 ); // retrieve the result from job #1
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that it wants the job with index 1 result #1", 1 ); // result index 1

		// finally, node 2 sends the job results from job 2
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"I have a message for you, archive\" B", 2 ); // node 2 says "I have a message for you, archive"
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"please archive this job result\" B", mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that the job id is 2", 2 ); // job_ib
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that the job id is 2 result index 1", 1 ); // result index
		std::string serialized_job_result2 = ts_assert_mpi_buffer_has_string( 2, "node 2 sends the serialized trpcage pose result for job 2" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_and_result2 = deserialized_larval_job_and_job_result( serialized_job_result2 );
		TS_ASSERT( job_and_result2.first );
		TS_ASSERT( job_and_result2.second );
		if ( job_and_result2.first ) {
			TS_ASSERT_EQUALS( job_and_result2.first->job_index(), 2 );
			TS_ASSERT_EQUALS( job_and_result2.first->nstruct_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->nstruct_max(), 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->input_job_result_indices().size(), 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->input_job_result_indices()[ 1 ].first, 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->input_job_result_indices()[ 1 ].second, 1 );
		}
		PoseJobResultOP pose_result2 = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_and_result2.second );
		TS_ASSERT( pose_result2 );
		if ( pose_result2 ) {
			TS_ASSERT_EQUALS( pose_result2->pose()->total_residue(), 20 );
		}

#endif
	}

	void test_MPIWorkPoolJobDistributor_master_node_vanilla_end_to_end_all_successes_jq_requires_larval_jobs_for_summaries() {
		TS_ASSERT( true ); // appease non-serialization builds


#ifdef SERIALIZATION

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 2 );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 1 starts its two part job-completion process,
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // result id
		send_string_to_node( 0, serialized_larval_job_and_job_result( 1, 1, 1, trpcage_pose_result ));

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_string_to_node( 0, serialized_larval_job( 1, 1, 1 ) ); // NEW! Node 1 has to send the serialized larval job back to node 0
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 1 asks for a new job
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// ok -- now we create a JQ on node 0, set it up to produce three jobs
		// and then create a job distributor for node 0 and tell it to go go go!

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ4OP jq( new PoolJQ4 );
		jq->njobs_ = 1;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int > ( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq->status_[ 1 ], jd3_job_status_success );

		// The job distributor should have used the note_job_completed function that takes
		// a LarvalJobOP as its first argument (and not the one that takes a Size as its first
		// argument
		TS_ASSERT_EQUALS( jq->jobs_completed_through_larval_job_interface_.size(), 1 );
		TS_ASSERT_EQUALS( jq->jobs_completed_through_larval_job_interface_[ 1 ], 1 );

		// The job distributor should have used the completed_job_summary function that takes
		// a LarvalJobOP as its first argument (and not the one that takes a Size as its first
		// argument
		TS_ASSERT_EQUALS( jq->summaries_through_larval_job_interface_.size(), 1 );
		TS_ASSERT_EQUALS( jq->summaries_through_larval_job_interface_[ 1 ], 1 );

		TS_ASSERT_EQUALS( jq->summaries_.count( sp1(1) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq->summaries_[ sp1(1) ])->energy(), 1234 );

		TS_ASSERT_EQUALS( jq->results_.count( sp1(1) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq->results_[ sp1(1) ])->pose()->total_residue(), 20 );

#endif
	}

	void test_workpool_jd3_worker_node_vanilla_end_to_end_all_successes_jq_requires_larval_jobs_for_summaries() {
		TS_ASSERT( true );


#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank(0);

		// send 1st job
		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// tell node 1 where to send the results of job 1:
		send_integers_to_node( 1, utility::vector1< int >( 1, 0 ) );
		send_integer_to_node( 1, mpi_work_pool_jd_archival_completed );

		// Now tell node 1 to spin dow
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		try {
			SimulateMPI::set_mpi_rank( 1 );
			PoolJQ4OP jq4( new PoolJQ4 );
			MPIWorkPoolJobDistributor jd;
			jd.go( jq4 );
		} catch ( ... ) {
			std::cerr << "Exception thrown from worker node in a unit test that should not have thrown an exception!" << std::endl;
			TS_ASSERT( false );
			return;
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Then it replies that its job is complete
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 has finished its job", mpi_work_pool_jd_job_success );
		ts_assert_mpi_buffer_has_size( 1, "node 1 finished job 1", 1 );
		ts_assert_mpi_buffer_has_size( 1, "node 1 finished job 1 producing 1 result", 1 );

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request C", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 has finished its job", mpi_work_pool_jd_archive_job_result );
		ts_assert_mpi_buffer_has_size( 1, "node 1 job 1", 1 );
		ts_assert_mpi_buffer_has_size( 1, "node 1 job 1 result 1", 1 );
		std::string serialized_job_result1 = ts_assert_mpi_buffer_has_string( 1, "node 1 send the serialized trpcage pose result" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_and_result1 = deserialized_larval_job_and_job_result( serialized_job_result1 );
		TS_ASSERT( job_and_result1.first );
		TS_ASSERT( job_and_result1.second );
		if ( job_and_result1.first ) {
			TS_ASSERT_EQUALS( job_and_result1.first->job_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_max(), 1 );
		}
		PoseJobResultOP pose_result1 = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_and_result1.second );
		TS_ASSERT( pose_result1 );
		if ( pose_result1 ) {
			TS_ASSERT_EQUALS( pose_result1->pose()->total_residue(), 20 );
		}

		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request D", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I finished archiving my completed job result -- I'm ready to send a JobSummary", mpi_work_pool_jd_job_success_and_archival_complete );
		std::string serialized_larval_job1 = ts_assert_mpi_buffer_has_string( 1, "node 1 finished this particular larval job (for job 1)" );
		LarvalJobOP larval_job1 = deserialize_larval_job( serialized_larval_job1 );
		TS_ASSERT( larval_job1 );
		if ( larval_job1 ) {
			TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
			TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
			TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 1 );
		}

		ts_assert_mpi_buffer_has_integer( 1, "node 1 is reporting the job status for job 1", jd3_job_status_success );
		std::string serialized_summary1 = ts_assert_mpi_buffer_has_string( 1, "node 1 says to node 0: here's the summary" );
		utility::vector1< JobSummaryOP > summaries1 = deserialized_job_summaries( serialized_summary1 );
		TS_ASSERT_EQUALS( summaries1.size(), 1 );
		if ( summaries1.size() != 1 ) return;
		TS_ASSERT( summaries1[1] );
		EnergyJobSummaryOP energy_summary1 = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( summaries1[1] );
		TS_ASSERT( energy_summary1 );
		if ( energy_summary1 ) {
			TS_ASSERT_EQUALS( energy_summary1->energy(), 1234 );
		}

		// OK -- now node 1 should ask for a job for a final time
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );


#endif
	}

	void test_master_node_response_to_failure_to_retrieve_job_result_from_archive()
	{
		TS_ASSERT( true );


#ifdef SERIALIZATION

		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		SimulateMPI::initialize_simulation( 3 );

		// ok -- the set up will be that node 1 is going to ask the archive for job result #1
		// and the archive is going to say that the job result isn't there.
		// In this test, the master node will think that the archive should have the job result.

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );


		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 2 starts its two part job-completion process,
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// and node 0 should say that node 2 ought to archive its result on node 1

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // node 1 says it finished job #1
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 2 asks for a new job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 2 requests a job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 1 tells node 0 about a failure to retrieve a job result
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_failed_to_retrieve_job_result );
		send_integer_to_node( 0, 2 ); // node 2 was the one who requested a job
		send_size_to_node( 0, 1 ); // it was job #1 that node 2 requested.
		send_size_to_node( 0, 1 ); // it was job #1 result #1 that node 2 requested.

		SimulateMPI::set_mpi_rank( 0 );

		PoolJQ3OP jq3( new PoolJQ3 );
		jq3->node_1_njobs_ = jq3->node_2_njobs_ = 1;
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3 );
			TS_ASSERT( false ); // this should not be reached.
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Failed to retrieve job result (1, 1) which was requested from node 1 by node 2"
				" but was not present there.\nJobDistributor on node 0 thinks the result should have been on node 1\n" );
		}

#endif
	}

	void test_master_node_response_to_failure_to_retrieve_job_result_from_archive_job_result_discarded()
	{
		TS_ASSERT( true );


#ifdef SERIALIZATION

		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		SimulateMPI::initialize_simulation( 4 );

		// there will be two rounds: round 1 has two jobs, round 2 has two jobs.
		// After round 1, job #1 will be preserved on the archive, and job #2 will be discarded.
		// Jobs 3 and 4 will list job #1's output as a required input.
		// Job 3 will run and finish, and then the job queen will tell the job distributor
		// to discard the result from job 1.
		// THEN Job 4 will request job 1 from the archive and the archive will tell node 0
		// that it doesn't have it.
		// The test is to make sure that the JobDistributor can provide a useful output
		// message saying that job 1 had been discarded.

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 2 starts its two part job-completion process,
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// and node 0 should say that node 2 ought to archive its result on node 1

		// node 3 starts its two part job-completion process,
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 2 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// and node 0 should say that node 2 ought to archive its result on node 1

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // node 2 says it finished job #1
		send_integer_to_node( 0, jd3_job_status_success ); // job #1 completed successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 2 ); // node 3 says it finished job #1
		send_integer_to_node( 0, jd3_job_status_success ); // job #2 completed successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 2 );
		// Now node 2 asks for a new job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		// Now node 3 asks for a new job
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );


		// now node 2 is going to say it has finished its job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 3 ); // node 1 says it finished job #3
		send_integer_to_node( 0, jd3_job_status_success ); // job #3 completed successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// And we pretend that node 3 is only just now getting started in retrieving
		// the job inputs from the archive but oh no! the job inputs it needs aren't there
		// any more
		// node 1 tells node 0 about a failure to retrieve a job result
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_failed_to_retrieve_job_result );
		send_integer_to_node( 0, 3 ); // node 3 was the one who requested a job
		send_size_to_node( 0, 1 ); // it was job #1 that node 3 requested.
		send_size_to_node( 0, 1 ); // it was job #1 result #1 that node 3 requested.

		SimulateMPI::set_mpi_rank( 0 );

		PoolJQ3bOP jq3b( new PoolJQ3b );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3b );
			TS_ASSERT( false ); // this should not be reached.
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Failed to retrieve job result (1, 1) which was requested from node 1 by node 3"
				" but was not present there.\nThis job is not listed as still running nor as having its JobResult stored on any archive node; it perhaps has been output or discarded already.\n" );
		}

#endif
	}

	void test_master_node_response_to_failure_to_retrieve_job_result_from_archive_job_previously_output()
	{
		TS_ASSERT( true );


#ifdef SERIALIZATION

		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		SimulateMPI::initialize_simulation( 4 );

		// there will be two rounds: round 1 has two jobs, round 2 has two jobs.
		// After round 1, job #1 will be preserved on the archive, and job #2 will be discarded.
		// Jobs 3 and 4 will list job #1's output as a required input.
		// Job 3 will run and finish, and then the job queen will tell the job distributor
		// to output the result from job 1 (removing it from the archive).
		// THEN Job 4 will request job 1 from the archive and the archive will tell node 0
		// that it doesn't have it.
		// The test is to make sure that the JobDistributor can provide a useful output
		// message saying that job 1 had been discarded.


		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 2 starts its two part job-completion process,
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// and node 0 should say that node 2 ought to archive its result on node 1

		// node 3 starts its two part job-completion process,
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 2 ); // job id
		send_size_to_node( 0, 1 ); // num results
		// and node 0 should say that node 2 ought to archive its result on node 1

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // node 2 says it finished job #1
		send_integer_to_node( 0, jd3_job_status_success ); // node 2 says it finished job #1 successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 2 ); // node 3 says it finished job #2
		send_integer_to_node( 0, jd3_job_status_success ); // node 3 says it finished job #2 successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 2 );
		// Now node 2 asks for a new job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		// Now node 3 asks for a new job
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );


		// now node 2 is going to say it has finished its job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 3 ); // node 1 says it finished job #3
		send_integer_to_node( 0, jd3_job_status_success ); // node 2 says it finished job #3 successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 1 );
		// ok, now the archive is going to receive a job result requests from
		// node 0 for job 1
		send_integer_to_node( 0, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( 0, serialized_larval_job_and_job_result( 2, 1, 1, trpcage_pose_result ));


		// And we pretend that node 3 is only just now getting started in retrieving
		// the job inputs from the archive .. so node 3 and node 1 start talking (but we don't get
		// to see their conversation) but oh no! the job inputs it needs aren't there on the archive
		// any more!
		// The archive then tells node 0 about a failure to retrieve a job result
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_failed_to_retrieve_job_result );
		send_integer_to_node( 0, 3 ); // node 3 was the one who requested a job
		send_size_to_node( 0, 1 ); // it was job #1 that node 3 requested.
		send_size_to_node( 0, 1 ); // it was job #1 result #1 that node 3 requested.

		SimulateMPI::set_mpi_rank( 0 );

		PoolJQ3bOP jq3b( new PoolJQ3b );
		jq3b->discard_job1_result_ = false; // instead, output job 1's result (removing it from the archive)
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3b );
			TS_ASSERT( false ); // this should not be reached.
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Failed to retrieve job result (1, 1) which was requested from node 1 by node 3"
				" but was not present there.\nThis job is not listed as still running nor as having its JobResult stored on any archive node; it perhaps has been output or discarded already.\n" );
		}

#endif
	}


	void test_work_pool_jd_archive_asked_for_result_it_doesnt_have()
	{
		TS_ASSERT( true );


#ifdef SERIALIZATION

		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		SimulateMPI::initialize_simulation( 4 );

		SimulateMPI::set_mpi_rank( 2 );

		// now retrieve the job result for job 1
		send_integer_to_node( 1, 2 );
		send_integer_to_node( 1, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 1, 1 ); // retrieve the result from job #1
		send_size_to_node( 1, 1 ); // retrieve the result from job #1

		// the archive doesn't directly exit -- it expects the JobDistributor on node 0 to exit once it has
		// recieved the error message (and when one MPI process exits, the whole job dies), so we actually
		// have to send a spin down signal to the archive to get the archive to exit it's listener loop
		SimulateMPI::set_mpi_rank( 0 );
		send_integer_to_node( 1, 0 );
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );


		// OK -- let's see what the archive node does
		SimulateMPI::set_mpi_rank( 1 );
		PoolJQ3bOP jq( new PoolJQ3b );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK -- what message did nodes 0 and 2 receive?
		SimulateMPI::set_mpi_rank( 0 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 tells node 0 that it needs to talk with it", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that it failed to retrieve a job result", mpi_work_pool_jd_failed_to_retrieve_job_result );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says that it was node 2 that asked for the job result", 2 ); // node 2 did the requesting
		ts_assert_mpi_buffer_has_size( 1, "node 1 says that node 2 was asking for job 1", 1 ); // it was job 1 whose result was requested
		ts_assert_mpi_buffer_has_size( 1, "node 1 says that node 2 was asking for job 1 result 1 ", 1 ); // it was job 1 result 1 that was requested

		SimulateMPI::set_mpi_rank( 2 );
		ts_assert_mpi_buffer_has_integer( 1, "node 2 says that it could not retrieve the result that node 2 had requested", mpi_work_pool_jd_failed_to_retrieve_job_result );

#endif
	}

	void test_work_pool_jd_manager_asked_for_result_it_doesnt_have()
	{
		TS_ASSERT( true );


#ifdef SERIALIZATION

		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 1 );

		// just cut right to the chase: before even running a single job, node 1 is going to ask node 0 for a job result
		// that node 0 does not have.
		// now retrieve the job result for job 1
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 0, 1 ); // retrieve the result from job #1
		send_size_to_node( 0, 1 ); // retrieve the result from job #1 result #1

		SimulateMPI::set_mpi_rank( 0 );

		PoolJQ3bOP jq3b( new PoolJQ3b );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3b );
			TS_ASSERT( false ); // this should not be reached.
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Failed to retrieve job result (1, 1) which was requested from node 0 by node 1"
				" but was not present there.\nThis job is not listed"
				" as still running nor as having its JobResult stored on any archive node; it perhaps has been output or discarded already.\n" );

		}
#endif
	}

	void test_send_bad_message_to_master_node_should_never_happen()
	{
		TS_ASSERT( true );


#ifdef SERIALIZATION

		SimulateMPI::initialize_simulation( 2 );

		SimulateMPI::set_mpi_rank( 1 );

		// just cut right to the chase: before even running a single job, node 1 is going to send a spin_down
		// signal to node 0, and node 0 just doesn't take that kind of instruction from a worker.
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 0 );

		PoolJQ3bOP jq3b( new PoolJQ3b );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq3b );
			TS_ASSERT( false ); // this should not be reached.
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Internal Error in the MPIWorkPoolJobDistributor: "
				"recieved inappropriate signal on node 0 from node 1: "
				+ utility::to_string( mpi_work_pool_jd_spin_down ) );
		}
#endif
	}

	void test_work_pool_jd_job_digraph_updater_and_job_node_w_zero_jobs()
	{
		TS_ASSERT( true );

#ifdef SERIALIZATION

		SimulateMPI::initialize_simulation( 2 );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 1 starts its two part job-completion process,
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // result id
		send_string_to_node( 0, serialized_larval_job_and_job_result( 1, 1, 1, trpcage_pose_result ));

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // job id
		send_integer_to_node( 0, jd3_job_status_success ); // node 1 says it finished job #1 successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 1 starts its two part job-completion process,
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 2 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 2 ); // job id
		send_size_to_node( 0, 1 ); // num results
		send_string_to_node( 0, serialized_larval_job_and_job_result( 1, 1, 2, trpcage_pose_result ));

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 2 ); // job id
		send_integer_to_node( 0, jd3_job_status_success ); // node 1 says it finished job #2 successfully
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 1 asks for a new job, but won't get one
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// Now let's create a JQ that will use the JobDigraphUpdater
		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ5OP jq( new PoolJQ5 );

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// And now we read the mail we got from node 0
		SimulateMPI::set_mpi_rank( 1 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv2 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job2 = deserialize_larval_job( ser_larv2 );
		TS_ASSERT_EQUALS( larval_job2->job_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_max(), 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );


#endif
	}

	void test_work_pool_jd_master_tries_to_serialize_larval_job_w_unserializable_data()
	{
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// Now let's create a JQ that creates unserializable larval jobs
		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ6OP jq( new PoolJQ6 );

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
			TS_ASSERT( false ); // we should not reach here
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::ostringstream oss;
			oss << "Failed to serialize LarvalJob 1 because it holds an unserializable object.  The following error message comes from the cereal library:\n" <<
				"Trying to save an unregistered polymorphic type (Unserializable).\n" <<
				"Make sure your type is registered with CEREAL_REGISTER_TYPE and that the archive you are" <<
				" using was included (and registered with CEREAL_REGISTER_ARCHIVE) prior to calling CEREAL_REGISTER_TYPE.\n" <<
				"If your type is already registered and you still see this error, you may need to use CEREAL_REGISTER_DYNAMIC_INIT.";

			//std::cout << "Caught exception\n" << e.msg() << "----\n";
			TS_ASSERT_EQUALS( e.msg(), oss.str() );
		}

#endif
	}

	void test_work_pool_jd_worker_tries_to_deserialize_larval_job_w_undeserializable_data()
	{
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		StandardInnerLarvalJobOP inner_job( new StandardInnerLarvalJob( 2 ) );
		inner_job->const_data_map().add( "testing", "testing", UndeserializableOP( new Undeserializable ));
		LarvalJobOP larval_job( new LarvalJob( inner_job, 1, 1 ));

		std::string undeserializable_larval_job;
		try {
			undeserializable_larval_job = serialized_T( larval_job );
		} catch ( cereal::Exception & e ) {
			TS_ASSERT( false );
			std::cerr << "Failed to serialize a LarvalJob holding a Undeserializable object\n" << e.what() << std::endl;
			return;
		}

// node 0 has a job for node 1
		SimulateMPI::set_mpi_rank( 0 );
		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, undeserializable_larval_job );

		// ok -- the worker node will be shut down when the mpi process exits; and this will happen when node 0
		// has output the error message that we sent and then exits; so for this unit test, just send node 1 a
		// spin-down signal
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		// Now let's create a JQ that creates unserializable larval jobs
		SimulateMPI::set_mpi_rank( 1 );
		PoolJQ1OP jq( new PoolJQ1 );

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false );
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Node 1 then says that it was unable to deserialize the larval job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_error );
		ts_assert_mpi_buffer_has_string( 1, "node 1 sent an error message",
			"Failed to deserialize larval job on worker node 1. Exiting.\nError message from cereal library:\n"
			"Undeserializable could not be deserialized\n" );

#endif
	}

	void test_work_pool_jd_worker_tries_to_deserialize_job_result_w_undeserializable_data()
	{
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		core::pose::PoseOP trpcage = create_trpcage_ideal_poseop();
		trpcage->set_const_data< Undeserializable >( "testing", "testing", Undeserializable());
		trpcage_pose_result->pose( trpcage );

		// node 0 has a job for node 1
		SimulateMPI::set_mpi_rank( 0 );
		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_T( create_larval_job_w_job1_dep( 1, 1, 2 ) ));
		send_sizes_to_node( 1, utility::vector1< Size >( 1, 0 ));

		// now node 0 sends an un-deserializable job result
		send_integer_to_node( 1, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( 1, serialized_larval_job_and_job_result( 1, 1, 1, trpcage_pose_result ));

		// ok -- the worker node will be shut down when the mpi process exits; and this will happen when node 0
		// has output the error message that we sent and then exits; so for this unit test, just send node 1 a
		// spin-down signal
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		// Now let's create a JQ that creates unserializable larval jobs
		SimulateMPI::set_mpi_rank( 1 );
		PoolJQ1OP jq( new PoolJQ1 );

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false );
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request B", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_retrieve_job_result );
		ts_assert_mpi_buffer_has_size( 1, "node 1 requests job #1", 1 );
		ts_assert_mpi_buffer_has_size( 1, "node 1 requests job #1 result #1", 1 );

		// Node 1 then says that it was unable to deserialize the job result
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_error );
		ts_assert_mpi_buffer_has_string( 1, "node 1 sent an error message",
			"Failed to deserialize LarvalJob & JobResult pair from job (1, 1) which is required as an input to job 2\n"
			"Error message from cereal library:\n"
			"Undeserializable could not be deserialized\n" );

#endif
	}

	void test_work_pool_jd_worker_tries_to_serialize_job_result_w_unserializable_data()
	{
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		// node 0 has a job for node 1
		SimulateMPI::set_mpi_rank( 0 );
		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// ok -- the worker node will be shut down when the mpi process exits; and this will happen when node 0
		// has output the error message that we sent and then exits; so for this unit test, just send node 1 a
		// spin-down signal
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		// Now let's create a JQ that creates unserializable larval jobs
		SimulateMPI::set_mpi_rank( 1 );
		PoolJQ8OP jq( new PoolJQ8 );

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false );
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Node 1 then says that it was unable to serialize the job result after job execution completed
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_error );
		ts_assert_mpi_buffer_has_string( 1, "node 1 sent an error message",
			"Failed to serialize LarvalJob or JobResult; all objects stored in these"
			" classes (or derived from them) must implement save & load serialization functions\nJob #1"
			" failed. Error message from the cereal library:\n"
			"Trying to save an unregistered polymorphic type (Unserializable).\n"
			"Make sure your type is registered with CEREAL_REGISTER_TYPE and that the archive you"
			" are using was included (and registered with CEREAL_REGISTER_ARCHIVE) prior to calling CEREAL_REGISTER_TYPE.\n"
			"If your type is already registered and you still see this error, you may need to use CEREAL_REGISTER_DYNAMIC_INIT.\n" );


#endif
	}

	void test_work_pool_jd_worker_tries_to_serialize_job_summary_w_unserializable_data()
	{
		TS_ASSERT( true );

#ifdef SERIALIZATION
		SimulateMPI::initialize_simulation( 2 );

		// node 0 has a job for node 1
		SimulateMPI::set_mpi_rank( 0 );
		send_integer_to_node( 1, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 1, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 1, utility::vector1< Size >() );

		// ok -- the worker node will be shut down when the mpi process exits; and this will happen when node 0
		// has output the error message that we sent and then exits; so for this unit test, just send node 1 a
		// spin-down signal
		send_integer_to_node( 1, mpi_work_pool_jd_spin_down );

		// Now let's create a JQ that creates unserializable larval jobs
		SimulateMPI::set_mpi_rank( 1 );
		PoolJQ9OP jq( new PoolJQ9 );

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
		} catch ( ... ) {
			TS_ASSERT( false );
		}

// NOW let's make sure that node 1 has been sending its messages the way it should have
		SimulateMPI::set_mpi_rank( 0 );
		// Node 1 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Node 1 then says that it was unable to deserialize the job summary
		ts_assert_mpi_buffer_has_integer( 1, "node 1 says to node 0: I have a request A", 1 );
		ts_assert_mpi_buffer_has_integer( 1, "node 1 requests a new job of node 0", mpi_work_pool_jd_error );
		ts_assert_mpi_buffer_has_string( 1, "node 1 sent an error message",
			"Failed to serialize JobSummary; all objects stored in this class (or that are derived from it) must implement save & load serialization functions\n"
			"Job #1 failed. Error message from the cereal library:\n"
			"Trying to save an unregistered polymorphic type (PoolJobSummary1).\n"
			"Make sure your type is registered with CEREAL_REGISTER_TYPE and that the archive you "
			"are using was included (and registered with CEREAL_REGISTER_ARCHIVE) prior to calling CEREAL_REGISTER_TYPE.\n"
			"If your type is already registered and you still see this error, you may need to use CEREAL_REGISTER_DYNAMIC_INIT.\n" );


#endif
	}

	void test_work_pool_jd_master_node_fails_to_deserialize_undeserializable_job_summary()
	{
		TS_ASSERT( true ); // appease non-serialization builds


#ifdef SERIALIZATION


		JobSummaryOP trpcage_job_summary( new PoolJobSummary2 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 2 );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 1 starts its two part job-completion process,
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 1, 1, 1, trpcage_pose_result ));

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// ok -- now we create a JQ on node 0, set it up to produce three jobs
		// and then create a job distributor for node 0 and tell it to go go go!

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ1OP jq( new PoolJQ1 );
		jq->njobs_ = 1;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cerr << e.msg() << std::endl;
			TS_ASSERT_EQUALS( e.msg(), "Failed to deserialize the JobSummary array for job #1\nError message from"
				" the cereal library:\nPoolJobSummary2 could not be deserialized\n" );
		}

#endif
	}

	void test_work_pool_jd_master_node_fails_to_deserialize_undeserializable_job_result()
	{
		TS_ASSERT( true ); // appease non-serialization builds

#ifdef SERIALIZATION

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );
		trpcage_pose_result->pose()->set_const_data< Undeserializable >( "testing", "testing", Undeserializable());

		SimulateMPI::initialize_simulation( 2 );

		// node 1 requests a job
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// node 1 starts its two part job-completion process,
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 1, 1, 1, trpcage_pose_result ));

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// ok -- now we create a JQ on node 0, set it up to produce three jobs
		// and then create a job distributor for node 0 and tell it to go go go!

		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ1OP jq( new PoolJQ1 );
		jq->njobs_ = 1;

		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			//std::cerr << e.msg() << std::endl;
			TS_ASSERT_EQUALS( e.msg(), "Failed to deserialize LarvalJob & JobResult pair for job #1 result index #1\nError message from"
				" the cereal library:\nUndeserializable could not be deserialized\n" );
		}

#endif
	}

	/////////////// Test the deallocation system
	void test_jd3_workpool_manager_deallocation_messages_sent_to_remote_nodes() {

#ifdef SERIALIZATION
		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 4 );

		// ok, all the nodes at once start requesting jobs from node 0
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );

		// let's start the job at node 1
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// and now the job at node 3
		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		// ok -- now let's pretend that node 2 has finished its job
		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );

		// ok -- before node 2 can finish its two part job-completion process,
		// node 1 will butt in and start its two part job-completion process,
		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 2 ); // job_id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 2 ); // job_id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 2, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 1 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 1, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 2 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 3 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 3, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 1 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 1 asks for a new job
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 2 );
		// Now node 2 asks for a new job
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );


		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 3 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		// Now node 3 asks for a new job
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 1 );

		// now retrieve the job result for job 1
		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 0, 1 ); // retrieve the result from job #1
		send_size_to_node( 0, 1 ); // retrieve the result from job #1 result #1

		SimulateMPI::set_mpi_rank( 2 );

		// now retrieve the job result for job 1
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 0, 1 ); // retrieve the result from job #1
		send_size_to_node( 0, 1 ); // retrieve the result from job #1 result #1

		SimulateMPI::set_mpi_rank( 3 );

		// now retrieve the job result for job 1
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_retrieve_job_result );
		send_size_to_node( 0, 1 ); // retrieve the result from job #1
		send_size_to_node( 0, 1 ); // retrieve the result from job #1 result #1

		// Now the second round of job results start trickling in

		SimulateMPI::set_mpi_rank( 1 );

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 4 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 4 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 1, 4, trpcage_pose_result ));

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 4 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 1 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 3 );
		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 6 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 6 ); // job id
		send_size_to_node( 0, 1 ); // num results
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 3, 6, trpcage_pose_result ));

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 6 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 3 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 2 );
		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success );
		send_size_to_node( 0, 5 ); // job id
		send_size_to_node( 0, 1 ); // num results

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_archive_job_result );
		send_size_to_node( 0, 5 ); // job id
		send_size_to_node( 0, 1 ); // result index
		send_string_to_node( 0, serialized_larval_job_and_job_result( 3, 2, 5, trpcage_pose_result ));

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_job_success_and_archival_complete );
		send_size_to_node( 0, 5 ); // job id
		send_integer_to_node( 0, jd3_job_status_success );
		send_string_to_node( 0, serialized_trpcage_job_summaries );

		send_integer_to_node( 0, 2 );
		send_integer_to_node( 0, mpi_work_pool_jd_new_job_request );

		//OK! Now create a PoolJQ11 and let 'er rip
		SimulateMPI::set_mpi_rank( 0 );
		PoolJQ11OP jq11( new PoolJQ11 );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq11 );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 0

		SimulateMPI::set_mpi_rank( 1 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv1 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job1 = deserialize_larval_job( ser_larv1 );
		TS_ASSERT_EQUALS( larval_job1->job_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job1->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are deallocations that should be made", mpi_work_pool_jd_deallocation_message );
		std::string node1_dealloc_msgs_str = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized deallocation message list" );
		std::list< deallocation::DeallocationMessageOP > node1_dealloc_msgs = deserialize_deallocation_msg_list( node1_dealloc_msgs_str );
		TS_ASSERT_EQUALS( node1_dealloc_msgs.size(), 2 );
		utility::vector1< deallocation::DeallocationMessageOP > node1_dealloc_msgs_vect( 2 );
		std::copy( node1_dealloc_msgs.begin(), node1_dealloc_msgs.end(), node1_dealloc_msgs_vect.begin() );
		deallocation::InputPoseDeallocationMessageOP n1msg1 = utility::pointer::dynamic_pointer_cast< deallocation::InputPoseDeallocationMessage > ( node1_dealloc_msgs_vect[1] );
		deallocation::InputPoseDeallocationMessageOP n1msg2 = utility::pointer::dynamic_pointer_cast< deallocation::InputPoseDeallocationMessage > ( node1_dealloc_msgs_vect[2] );
		TS_ASSERT( n1msg1 );
		TS_ASSERT( n1msg2 );
		TS_ASSERT_EQUALS( n1msg1->pose_id(), 1 );
		TS_ASSERT_EQUALS( n1msg2->pose_id(), 2 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv4 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 1 the serialized LarvalJob" );
		LarvalJobOP larval_job4 = deserialize_larval_job( ser_larv4 );
		TS_ASSERT_EQUALS( larval_job4->job_index(), 4 );
		TS_ASSERT_EQUALS( larval_job4->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( larval_job4->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job4->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 1 a vector1 of job result indices with one element whose value is 0", utility::vector1< Size >( 1, 0 ) );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1a = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1a = deserialized_larval_job_and_job_result( ser_job_res1a );
		TS_ASSERT( job_res1a.first );
		TS_ASSERT( job_res1a.second );
		TS_ASSERT_EQUALS( job_res1a.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1a.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1a.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1a.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 1 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 1 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 2 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv2 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized LarvalJob" );
		LarvalJobOP larval_job2 = deserialize_larval_job( ser_larv2 );
		TS_ASSERT_EQUALS( larval_job2->job_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( larval_job2->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 2 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 2 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that there are deallocations that should be made", mpi_work_pool_jd_deallocation_message );
		std::string node2_dealloc_msgs_str = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized deallocation message list" );
		std::list< deallocation::DeallocationMessageOP > node2_dealloc_msgs = deserialize_deallocation_msg_list( node2_dealloc_msgs_str );
		TS_ASSERT_EQUALS( node2_dealloc_msgs.size(), 2 );
		utility::vector1< deallocation::DeallocationMessageOP > node2_dealloc_msgs_vect( 2 );
		std::copy( node2_dealloc_msgs.begin(), node2_dealloc_msgs.end(), node2_dealloc_msgs_vect.begin() );
		deallocation::InputPoseDeallocationMessageOP n2msg1 = utility::pointer::dynamic_pointer_cast< deallocation::InputPoseDeallocationMessage > ( node2_dealloc_msgs_vect[1] );
		deallocation::InputPoseDeallocationMessageOP n2msg2 = utility::pointer::dynamic_pointer_cast< deallocation::InputPoseDeallocationMessage > ( node2_dealloc_msgs_vect[2] );
		TS_ASSERT( n2msg1 );
		TS_ASSERT( n2msg2 );
		TS_ASSERT_EQUALS( n2msg1->pose_id(), 1 );
		TS_ASSERT_EQUALS( n2msg2->pose_id(), 2 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv5 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 2 the serialized LarvalJob" );
		LarvalJobOP larval_job5 = deserialize_larval_job( ser_larv5 );
		TS_ASSERT_EQUALS( larval_job5->job_index(), 5 );
		TS_ASSERT_EQUALS( larval_job5->nstruct_index(), 2 );
		TS_ASSERT_EQUALS( larval_job5->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job5->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 2 a vector1 of job result indices with one element whose value is 0", utility::vector1< Size >( 1, 0 ) );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1b = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1b = deserialized_larval_job_and_job_result( ser_job_res1b );
		TS_ASSERT( job_res1b.first );
		TS_ASSERT( job_res1b.second );
		TS_ASSERT_EQUALS( job_res1b.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1b.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1b.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1b.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 2 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 2 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		SimulateMPI::set_mpi_rank( 3 );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has a job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv3 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized LarvalJob" );
		LarvalJobOP larval_job3 = deserialize_larval_job( ser_larv3 );
		TS_ASSERT_EQUALS( larval_job3->job_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( larval_job3->nstruct_max(), 3 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 3 an empty list of job result indices", utility::vector1< Size >() );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 3 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that there are deallocations that should be made", mpi_work_pool_jd_deallocation_message );
		std::string node3_dealloc_msgs_str = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized deallocation message list" );
		std::list< deallocation::DeallocationMessageOP > node3_dealloc_msgs = deserialize_deallocation_msg_list( node3_dealloc_msgs_str );
		TS_ASSERT_EQUALS( node3_dealloc_msgs.size(), 2 );
		utility::vector1< deallocation::DeallocationMessageOP > node3_dealloc_msgs_vect( 2 );
		std::copy( node3_dealloc_msgs.begin(), node3_dealloc_msgs.end(), node3_dealloc_msgs_vect.begin() );
		deallocation::InputPoseDeallocationMessageOP n3msg1 = utility::pointer::dynamic_pointer_cast< deallocation::InputPoseDeallocationMessage > ( node3_dealloc_msgs_vect[1] );
		deallocation::InputPoseDeallocationMessageOP n3msg2 = utility::pointer::dynamic_pointer_cast< deallocation::InputPoseDeallocationMessage > ( node3_dealloc_msgs_vect[2] );
		TS_ASSERT( n3msg1 );
		TS_ASSERT( n3msg2 );
		TS_ASSERT_EQUALS( n3msg1->pose_id(), 1 );
		TS_ASSERT_EQUALS( n3msg2->pose_id(), 2 );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has another job for it to run", mpi_work_pool_jd_new_job_available );
		std::string ser_larv6 = ts_assert_mpi_buffer_has_string( 0, "node 0 sends node 3 the serialized LarvalJob" );
		LarvalJobOP larval_job6 = deserialize_larval_job( ser_larv6 );
		TS_ASSERT_EQUALS( larval_job6->job_index(), 6 );
		TS_ASSERT_EQUALS( larval_job6->nstruct_index(), 3 );
		TS_ASSERT_EQUALS( larval_job6->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices().size(), 1 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices()[ 1 ].first, 1 );
		TS_ASSERT_EQUALS( larval_job6->input_job_result_indices()[ 1 ].second, 1 );
		ts_assert_mpi_buffer_has_sizes( 0, "node 0 sends node 3 a vector1 of job result indices with one element whose value is 0", utility::vector1< Size >( 1, 0 ) );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that it has the job result it was expecting", mpi_work_pool_jd_job_result_retrieved );
		std::string ser_job_res1c = ts_assert_mpi_buffer_has_string( 0, "node 0 sends the serialized larval-job and job-result from job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_res1c = deserialized_larval_job_and_job_result( ser_job_res1b );
		TS_ASSERT( job_res1c.first );
		TS_ASSERT( job_res1c.second );
		TS_ASSERT_EQUALS( job_res1c.first->job_index(), 1 );
		TS_ASSERT_EQUALS( job_res1c.first->nstruct_index(), 1 );
		TS_ASSERT_EQUALS( job_res1c.first->nstruct_max(), 3 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_res1c.second )->pose()->total_residue(), 20 );

		ts_assert_mpi_buffer_has_integers( 0, "node 0 tells node 3 to archive its results on node 0", utility::vector1< int >( 1, 0 ) );
		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that archival completed", mpi_work_pool_jd_archival_completed );

		ts_assert_mpi_buffer_has_integer( 0, "node 0 tells node 3 that there are no jobs left to run", mpi_work_pool_jd_spin_down );

		TS_ASSERT_EQUALS( jq11->status_.count( 1 ), 1 );
		TS_ASSERT_EQUALS( jq11->status_.count( 2 ), 1 );
		TS_ASSERT_EQUALS( jq11->status_.count( 3 ), 1 );
		TS_ASSERT_EQUALS( jq11->status_.count( 4 ), 1 );
		TS_ASSERT_EQUALS( jq11->status_.count( 5 ), 1 );
		TS_ASSERT_EQUALS( jq11->status_.count( 6 ), 1 );
		TS_ASSERT_EQUALS( jq11->status_[ 1 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq11->status_[ 2 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq11->status_[ 3 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq11->status_[ 4 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq11->status_[ 5 ], jd3_job_status_success );
		TS_ASSERT_EQUALS( jq11->status_[ 6 ], jd3_job_status_success );

		TS_ASSERT_EQUALS( jq11->summaries_.count( sp1( 1 ) ), 1 );
		TS_ASSERT_EQUALS( jq11->summaries_.count( sp1( 2 ) ), 1 );
		TS_ASSERT_EQUALS( jq11->summaries_.count( sp1( 3 ) ), 1 );
		TS_ASSERT_EQUALS( jq11->summaries_.count( sp1( 4 ) ), 1 );
		TS_ASSERT_EQUALS( jq11->summaries_.count( sp1( 5 ) ), 1 );
		TS_ASSERT_EQUALS( jq11->summaries_.count( sp1( 6 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq11->summaries_[ sp1( 1 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq11->summaries_[ sp1( 2 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq11->summaries_[ sp1( 3 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq11->summaries_[ sp1( 4 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq11->summaries_[ sp1( 5 ) ])->energy(), 1234 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( jq11->summaries_[ sp1( 6 ) ])->energy(), 1234 );

		TS_ASSERT_EQUALS( jq11->results_.count( sp1( 1 ) ), 0 );
		TS_ASSERT_EQUALS( jq11->results_.count( sp1( 2 ) ), 0 );
		TS_ASSERT_EQUALS( jq11->results_.count( sp1( 3 ) ), 0 );
		TS_ASSERT_EQUALS( jq11->results_.count( sp1( 4 ) ), 1 );
		TS_ASSERT_EQUALS( jq11->results_.count( sp1( 5 ) ), 1 );
		TS_ASSERT_EQUALS( jq11->results_.count( sp1( 6 ) ), 1 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq11->results_[ sp1( 4 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq11->results_[ sp1( 5 ) ])->pose()->total_residue(), 20 );
		TS_ASSERT_EQUALS( utility::pointer::dynamic_pointer_cast< PoseJobResult > ( jq11->results_[ sp1( 6 ) ])->pose()->total_residue(), 20 );



#endif
	}

	void test_jd3_workpool_worker_deallocation_messages_sent_to_remote_nodes() {
		TS_ASSERT( true );

#ifdef SERIALIZATION
		core_init_with_additional_options( "-jd3::n_archive_nodes 1 -jd3::do_not_archive_on_node0 -jd3::compress_job_results 0" );

		EnergyJobSummaryOP trpcage_job_summary( new EnergyJobSummary );
		trpcage_job_summary->energy( 1234 );
		utility::vector1< JobSummaryOP > summaries( 1, trpcage_job_summary );
		std::string serialized_trpcage_job_summaries = serialized_job_summaries( summaries );

		PoseJobResultOP trpcage_pose_result( new PoseJobResult );
		trpcage_pose_result->pose(  create_trpcage_ideal_poseop() );

		SimulateMPI::initialize_simulation( 3 );

		SimulateMPI::set_mpi_rank( 0 );

		// send 1st job
		send_integer_to_node( 2, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 2, serialized_larval_job( 1, 1, 1 ) );
		send_sizes_to_node( 2, utility::vector1< Size >() );

		// tell node 2 where to send the results of job 1:
		send_integers_to_node( 2, utility::vector1< int >( 1, 1 ) );

		SimulateMPI::set_mpi_rank( 1 );
		// ok, node 2 will say that it needs to send us a job, and it'll send it to node 1 (the archive)
		// and then the archive needs to reply that the archival was successful
		send_integer_to_node( 2, mpi_work_pool_jd_archival_completed );

		SimulateMPI::set_mpi_rank( 0 );
		// send deallocation message
		send_integer_to_node( 2, mpi_work_pool_jd_deallocation_message );
		std::list< core::Size > pose_ids{ 1, 2 };
		send_string_to_node( 2, serialized_input_pose_deallocation_msg_list( pose_ids ) );

		// send 2nd job, which lists job 1 as a required input
		send_integer_to_node( 2, mpi_work_pool_jd_new_job_available );
		send_string_to_node( 2, serialized_T( create_larval_job_w_job1_dep( 1, 1, 2 ) ));
		send_sizes_to_node( 2, utility::vector1< Size >( 1, 1 ) ); // say that job 1's results live on node 1

		SimulateMPI::set_mpi_rank( 1 );
		send_integer_to_node( 2, mpi_work_pool_jd_job_result_retrieved );
		send_string_to_node( 2, serialized_larval_job_and_job_result( 1, 1, 1, trpcage_pose_result ));

		SimulateMPI::set_mpi_rank( 0 );
		// tell node 2 where to send the results of job 2:
		send_integers_to_node( 2, utility::vector1< int >( 1, 1 ) );

		SimulateMPI::set_mpi_rank( 1 );
		// ok, node 2 will say that it needs to send us a job, and it'll send it to node 1 (the archive)
		// and then the archive needs to reply that the archival was successful
		send_integer_to_node( 2, mpi_work_pool_jd_archival_completed );

		SimulateMPI::set_mpi_rank( 0 );
		send_integer_to_node( 2, mpi_work_pool_jd_spin_down ); // no more jobs to hand out

		//OK! Now create a PoolJQ11 and let 'er rip
		SimulateMPI::set_mpi_rank( 2 );
		PoolJQ11OP jq11( new PoolJQ11 );
		MPIWorkPoolJobDistributor jd;
		try {
			jd.go( jq11 );
		} catch ( ... ) {
			TS_ASSERT( false /*no exception should be thrown*/ );
		}

// OK!
// Now we read the mail we got from node 2

		SimulateMPI::set_mpi_rank( 0 );
		// Node 2 starts by requesting a job
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request A", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Then it replies that its job is complete
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request B", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 has finished its job", mpi_work_pool_jd_job_success );
		ts_assert_mpi_buffer_has_size( 2, "node 2 finished job #1", 1 );
		ts_assert_mpi_buffer_has_size( 2, "node 2 finished job #1 with 1 result", 1 );

		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request C", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I finished archiving my completed job result -- I'm ready to send a JobSummary", mpi_work_pool_jd_job_success_and_archival_complete );
		ts_assert_mpi_buffer_has_size( 2, "node 2 says to node 0: this was for job index 1", 1 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: job 1 was successful", jd3_job_status_success );
		std::string serialized_summary1 = ts_assert_mpi_buffer_has_string( 2, "node 2 says to node 0: here's the summary" );
		utility::vector1< JobSummaryOP > summaries1 = deserialized_job_summaries( serialized_summary1 );
		TS_ASSERT_EQUALS( summaries1.size(), 1 );
		if ( summaries1.size() != 1 ) return;
		TS_ASSERT( summaries1[ 1 ] );
		EnergyJobSummaryOP energy_summary1 = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( summaries1[ 1 ] );
		TS_ASSERT( energy_summary1 );
		if ( energy_summary1 ) {
			TS_ASSERT_EQUALS( energy_summary1->energy(), 1234 );
		}

		// Now Node 2 requests another job
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request D", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		// Then it replies that its job is complete
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request E", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 has finished its job", mpi_work_pool_jd_job_success );
		ts_assert_mpi_buffer_has_size( 2, "node 2 finished job #2", 2 );
		ts_assert_mpi_buffer_has_size( 2, "node 2 finished job #2 which generated one result", 1 );

		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request F", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I finished archiving my completed job result -- I'm ready to send a JobSummary", mpi_work_pool_jd_job_success_and_archival_complete );
		ts_assert_mpi_buffer_has_size( 2, "node 2 says to node 0: this was for job index 2", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: job 2 was successful", jd3_job_status_success );
		std::string serialized_summary2 = ts_assert_mpi_buffer_has_string( 2, "node 2 says to node 0: here's the summary" );
		utility::vector1< JobSummaryOP > summaries2 = deserialized_job_summaries( serialized_summary2 );
		TS_ASSERT_EQUALS( summaries2.size(), 1 );
		if ( summaries2.size() != 1 ) return;
		TS_ASSERT( summaries2[1] );
		EnergyJobSummaryOP energy_summary2 = utility::pointer::dynamic_pointer_cast< EnergyJobSummary > ( summaries2[1] );
		TS_ASSERT( energy_summary2 );
		if ( energy_summary2 ) {
			TS_ASSERT_EQUALS( energy_summary2->energy(), 1234 );
		}

		// Now Node 2 requests another job, and it will turn out that there are no more, so node 0 will send a spin down signal.
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to node 0: I have a request G", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 requests a new job of node 0", mpi_work_pool_jd_new_job_request );

		SimulateMPI::set_mpi_rank( 1 );

		// Let's see what node 2 sent to the archive node.
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"I have a message for you, archive\"", 2 ); // node 2 says "I have a message for you, archive"
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"please archive this job result\"", mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		ts_assert_mpi_buffer_has_size( 2, "node 2 says this job result is for job #1", 1 ); // job_id
		ts_assert_mpi_buffer_has_size( 2, "node 2 says this job result is for job #1 result #1", 1 ); // result index
		std::string serialized_job_result1 = ts_assert_mpi_buffer_has_string( 2, "node 2 sends the serialized trpcage pose result for job 1" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_and_result1 = deserialized_larval_job_and_job_result( serialized_job_result1 );
		TS_ASSERT( job_and_result1.first );
		TS_ASSERT( job_and_result1.second );
		if ( job_and_result1.first ) {
			TS_ASSERT_EQUALS( job_and_result1.first->job_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result1.first->nstruct_max(), 1 );
		}
		PoseJobResultOP pose_result1 = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_and_result1.second );
		TS_ASSERT( pose_result1 );
		if ( pose_result1 ) {
			TS_ASSERT_EQUALS( pose_result1->pose()->total_residue(), 20 );
		}

		// node 2 then asks for the job results from job 1
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says to the archive that it has a request", 2 );
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says that it needs a job result",  mpi_work_pool_jd_retrieve_job_result );
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that it wants the job 1", 1 ); // retrieve the result from job #1
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that it wants the job 1 with result index 1", 1 ); // retrieve the result from job #1 result #1

		// finally, node 2 sends the job results from job 2
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"I have a message for you, archive\" B", 2 ); // node 2 says "I have a message for you, archive"
		ts_assert_mpi_buffer_has_integer( 2, "node 2 says \"please archive this job result\" B", mpi_work_pool_jd_archive_job_result ); // "please archive this job result"
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that the job id is 2", 2 ); // job_ib
		ts_assert_mpi_buffer_has_size( 2, "node 2 says that the job id is 2", 1 ); // result index 1
		std::string serialized_job_result2 = ts_assert_mpi_buffer_has_string( 2, "node 2 sends the serialized trpcage pose result for job 2" );
		MPIWorkPoolJobDistributor::LarvalJobAndResult job_and_result2 = deserialized_larval_job_and_job_result( serialized_job_result2 );
		TS_ASSERT( job_and_result2.first );
		TS_ASSERT( job_and_result2.second );
		if ( job_and_result2.first ) {
			TS_ASSERT_EQUALS( job_and_result2.first->job_index(), 2 );
			TS_ASSERT_EQUALS( job_and_result2.first->nstruct_index(), 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->nstruct_max(), 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->input_job_result_indices().size(), 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->input_job_result_indices()[ 1 ].first, 1 );
			TS_ASSERT_EQUALS( job_and_result2.first->input_job_result_indices()[ 1 ].second, 1 );
		}
		PoseJobResultOP pose_result2 = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( job_and_result2.second );
		TS_ASSERT( pose_result2 );
		if ( pose_result2 ) {
			TS_ASSERT_EQUALS( pose_result2->pose()->total_residue(), 20 );
		}

		TS_ASSERT_EQUALS( jq11->received_messages_.size(), 2 );
		typedef deallocation::InputPoseDeallocationMessageOP DeallocPoseOP;
		typedef deallocation::InputPoseDeallocationMessage DeallocPose;
		DeallocPoseOP msg1 = utility::pointer::dynamic_pointer_cast< DeallocPose >( jq11->received_messages_[ 1 ] );
		DeallocPoseOP msg2 = utility::pointer::dynamic_pointer_cast< DeallocPose >( jq11->received_messages_[ 2 ] );
		TS_ASSERT( msg1 );
		TS_ASSERT( msg2 );
		TS_ASSERT_EQUALS( msg1->pose_id(), 1 );
		TS_ASSERT_EQUALS( msg2->pose_id(), 2 );

#endif
	}


};
