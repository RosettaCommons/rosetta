// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/MultiThreadedJobDistributor.cc
/// @brief  MultiThreadedJobDistributor class method definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef MULTI_THREADED

// Unit headers
#include <protocols/jd3/job_distributors/MultiThreadedJobDistributor.hh>

// Package headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/job_distributors/JobExtractor.hh>
#include <protocols/jd3/output/OutputSpecification.hh>
#include <protocols/jd3/output/ResultOutputter.hh>

// Core headers
#include <core/init/init.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <chrono>
#include <ctime>
#include <functional>
#include <string>
#include <thread>

// CTPL
#include <CTPL/ctpl_stl.h>

namespace protocols {
namespace jd3 {
namespace job_distributors {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.MultiThreadedJobDistributor" );

MultiThreadedJobDistributor::MultiThreadedJobDistributor() :
	job_extractor_( new JobExtractor ),
	nthreads_( 0 ),
	default_retry_limit_( 1 )
{
	if ( ! basic::options::option[ basic::options::OptionKeys::jd3::nthreads ].active() ) {
		throw utility::excn::EXCN_Msg_Exception( "jd3::nthreads not specified" );
	}

	if ( basic::options::option[ basic::options::OptionKeys::jd3::nthreads ]() <= 0 ) {
		throw utility::excn::EXCN_Msg_Exception( "negative or zero value passed in for jd3::nthreads" );
	}

	nthreads_ = basic::options::option[ basic::options::OptionKeys::jd3::nthreads ]();
	thread_pool_.reset( new ctpl::thread_pool( nthreads_ ) );
}

MultiThreadedJobDistributor::~MultiThreadedJobDistributor() {
	thread_pool_->stop( true );
}

void
MultiThreadedJobDistributor::go( JobQueenOP queen )
{
	job_queen_ = queen;
	job_extractor_->set_job_queen( queen );
	job_dag_ = job_extractor_->create_initial_job_dag();

	// extract the first nthreads_ jobs from the JobQueen, and start them running.
	core::Size n_queued( 0 );
	while ( n_queued < nthreads_ ) {
		if ( job_extractor_->job_queue_empty() ) break;
		LarvalJobOP larval_job = job_extractor_->pop_job_from_queue();
		//TR << "Popping job " << larval_job->job_index() << std::endl;
		if ( job_queen_->has_job_completed( larval_job ) ) {
			job_queen_->note_job_completed( larval_job, jd3_job_previously_executed, 0 );
			job_extractor_->note_job_no_longer_running( larval_job->job_index() );
			continue;
		}
		//TR << "Preparing next job: " << larval_job->job_index() << std::endl;
		if ( prepare_next_job( larval_job, 1 ) ) ++n_queued;
	}

	// Then ask the job_extractor_ for jobs as long as it says it is not_done()
	// and also scan over the list of running jobs to remove the ones that have completed,
	// giving the job results to the job queen
	while ( job_extractor_->not_done() ) {

		// Queue more work if needed
		while ( jobs_running_.size() <= 2*nthreads_ ) {
			LarvalJobOP larval_job;
			while ( ! job_extractor_->job_queue_empty() ) {
				larval_job = job_extractor_->pop_job_from_queue();
				//TR << "Popping job " << larval_job->job_index() << std::endl;
				if ( job_queen_->has_job_completed( larval_job ) ) {
					//TR << "Job " << larval_job->job_index() << " has already completed" << std::endl;
					job_queen_->note_job_completed( larval_job, jd3_job_previously_executed, 0 );
					job_extractor_->note_job_no_longer_running( larval_job->job_index() );
					larval_job = nullptr;
					continue;
				}
				break; // we have a job, so let's submit it
			}
			if ( ! larval_job ) break;
			//TR << "Preparing next job: " << larval_job->job_index() << std::endl;
			prepare_next_job( larval_job, 1 );
		}

		// Scan the running jobs -- process and delete completed elements from the list
		// as they are encountered
		bool job_finished( false );
		for ( auto job_runner_iter = jobs_running_.begin(); job_runner_iter != jobs_running_.end(); /* empty! */ ) {
			JobRunnerOP job_runner = *job_runner_iter;
			if ( job_runner->complete() ) {
				process_completed_job( job_runner );
				auto to_delete = job_runner_iter;
				++job_runner_iter;
				jobs_running_.erase( to_delete );
				job_finished = true;
			} else {
				++job_runner_iter;
			}
		}
		if ( ! job_extractor_->job_queue_empty() && ! job_finished ) {
			std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );
		}
	}

}


bool
MultiThreadedJobDistributor::prepare_next_job( LarvalJobOP larval_job, core::Size attempt_count )
{
	utility::vector1< JobResultCOP > input_job_results( larval_job->input_job_result_indices().size() );
	for ( Size ii = 1; ii <= input_job_results.size(); ++ii ) {
		auto iter = job_results_.find( larval_job->input_job_result_indices()[ ii ] );
		if ( iter == job_results_.end() ) {
			std::ostringstream oss;
			oss << "Failed to retrieve requested input for job #" << larval_job->job_index() << " of (";
			oss << larval_job->input_job_result_indices()[ ii ].first << ", ";
			oss << larval_job->input_job_result_indices()[ ii ].second << ")";
			TR.Error << oss.str() << std::endl;
			throw utility::excn::EXCN_Msg_Exception( oss.str() );
		}
		input_job_results[ ii ] = iter->second.second;
	}

	JobOP mature_job;
	try {
		mature_job = job_queen_->mature_larval_job( larval_job, input_job_results );
	} catch ( ... ) {
		job_queen_->note_job_completed( larval_job, jd3_job_status_inputs_were_bad, 0 );
		job_extractor_->note_job_no_longer_running( larval_job->job_index() );
		return false;
	}
	Size num_retries = larval_job->defines_retry_limit() ? larval_job->retry_limit() : default_retry_limit_;
	JobRunnerOP runner( new JobRunner( larval_job, mature_job, attempt_count, num_retries ));
	thread_pool_->push( std::bind( &JobRunner::run, runner, std::placeholders::_1 ));
	jobs_running_.push_back( runner );
	return true;
}


void
MultiThreadedJobDistributor::process_completed_job( JobRunnerOP job_runner )
{
	TR << "Finished job " << job_runner->larval_job()->job_index() << " in thread " << job_runner->running_thread() << std::endl;

	LarvalJobOP larval_job = job_runner->larval_job();
	CompletedJobOutput const & output = job_runner->job_output();
	// now we process the job
	if ( job_runner->exited_w_exception() ) {
		TR.Error << "Job exited with exception\n" << job_runner->exception_message() << std::endl;
		job_queen_->note_job_completed( larval_job, jd3_job_status_failed_w_exception, 0 );
	} else if ( output.status == jd3_job_status_failed_retry ) {
		Size num_retries = larval_job->defines_retry_limit() ? larval_job->retry_limit() : default_retry_limit_;
		if ( job_runner->attempt_count() >= num_retries ) {
			job_queen_->note_job_completed( larval_job, jd3_job_status_failed_max_retries, 0 );
		} else {
			// put the job back into the queue to try again
			prepare_next_job( larval_job, job_runner->attempt_count() + 1 );
			return;
		}
	} else {
		Size job_index = larval_job->job_index();
		job_queen_->note_job_completed( larval_job, output.status, output.job_results.size() );
		for ( Size ii = 1; ii <= output.job_results.size(); ++ii ) {
			job_queen_->completed_job_summary( larval_job, ii, output.job_results[ ii ].first );
			job_results_[ { job_index, ii } ] = std::make_pair( larval_job, output.job_results[ ii ].second );
		}
	}

	job_extractor_->note_job_no_longer_running( larval_job->job_index() );

	potentially_output_some_job_results();
	potentially_discard_some_job_results();

}

void
MultiThreadedJobDistributor::potentially_output_some_job_results()
{

	std::list< output::OutputSpecificationOP > jobs_to_output =
		job_queen_->jobs_that_should_be_output();
	for ( auto const output_spec : jobs_to_output ) {
		JobResultID result_id = output_spec->result_id();
		auto result_iter = job_results_.find( result_id );
		if ( result_iter == job_results_.end() ) {
			throw utility::excn::EXCN_Msg_Exception( "Failed to retrieve job result (" +
				utility::to_string( result_id.first )  + ", " + utility::to_string( result_id.second ) +
				") for outputting as requested by the JobQeen. Has this job already been output?" );
		}
		std::pair< LarvalJobOP, JobResultOP > const & job_to_output = result_iter->second;
		output::ResultOutputterOP outputter = job_queen_->result_outputter( *output_spec );
		outputter->write_output( *output_spec, *job_to_output.second );

		// and now erase the job
		job_results_.erase( result_id );
	}
}

void
MultiThreadedJobDistributor::potentially_discard_some_job_results()
{
	// ask the job queen if she wants to discard any results
	std::list< JobResultID >  jobs_to_discard = job_queen_->job_results_that_should_be_discarded();
	if ( jobs_to_discard.empty() ) return;

	for ( JobResultID result_id : jobs_to_discard ) {
		if ( job_results_.count( result_id ) == 0 ) {
			// TO DO: better diagnostic message here
			throw utility::excn::EXCN_Msg_Exception( "Failed to find job result " +
				utility::to_string( result_id.first ) + ", " + utility::to_string( result_id.second ) +
				" for discarding as requested by the JobQeen. Has this job already been" +
				" output or discarded?" );
		}

		job_results_.erase( result_id );
	}
}


JobRunner::JobRunner( LarvalJobOP larval_job, JobOP mature_job, core::Size attempt_count, core::Size retry_limit ) :
	attempt_count_( attempt_count ),
	retry_limit_( retry_limit ),
	larval_job_( larval_job ),
	mature_job_( mature_job ),
	exited_w_exception_( false ),
	complete_( false )
{}

void JobRunner::run( int thread_index )
{
	running_thread_ = thread_index;
	// The RNG needs to be set up for the thread

	core::init::RandomGeneratorSettings rgs;
	rgs.initialize_from_options( basic::options::option );

	if ( ! rgs.const_seed() ) {
		if ( ! numeric::random::rg().initialized() ) {
			rgs.seed_offset( rgs.seed_offset() + thread_index );
			rgs.mpi_bcast_seed_from_node0( false );
			int seed = core::init::determine_random_number_seed( rgs );
			core::init::init_random_generators( seed, rgs.rng_type() );
		}
	} else {
		core::init::init_random_generators(
			rgs.seed() + rgs.seed_offset() + retry_limit_ * larval_job_->job_index() + attempt_count_,
			rgs.rng_type() );
	}

	try {
		job_result_ = mature_job_->run();
	} catch ( utility::excn::EXCN_Base & e ) {
		exception_message_ = e.msg();
		exited_w_exception_ = true;
	} catch ( ... ) {
		exited_w_exception_ = true;
	}
	complete_.store( true );
}

bool JobRunner::exited_w_exception() const
{
	debug_assert( complete_.load() );
	return exited_w_exception_;
}

std::string const &
JobRunner::exception_message() const
{
	debug_assert( complete_.load() );
	return exception_message_;
}

LarvalJobOP JobRunner::larval_job() const { return larval_job_; }
JobOP JobRunner::mature_job() const { return mature_job_; }
core::Size JobRunner::attempt_count() const { return attempt_count_; }


CompletedJobOutput const &
JobRunner::job_output() const { return job_result_; }

bool JobRunner::complete() const { return complete_.load(); }

int JobRunner::running_thread() const { return running_thread_; }

} // job_distributors
} // jd3
} // protocols

#endif
