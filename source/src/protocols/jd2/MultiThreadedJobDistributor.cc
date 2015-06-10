// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MultiThreadedJobDistributor.cc
/// @brief  Job distributor that launches threads to carry out independant jobs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifdef MULTI_THREADED
#ifdef CXX11

// Unit headers
#include <protocols/jd2/MultiThreadedJobDistributor.hh>

//Package headers
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/Parser.hh>


// Project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/exit.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <utility/assert.hh>
#include <string>


// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <utility/vector1.hh>

// Boost includes
#include <boost/bind.hpp>


// C++11 headers
#include <chrono>
#include <thread>

using basic::Warning;
static thread_local basic::Tracer TR( "protocols.jd2.MultiThreadedJobDistributor" );

namespace protocols {
namespace jd2 {


thread_local RunningJobOP running_job_;

void
set_running_job( RunningJob const & job ) {
	running_job_ = RunningJobOP( new RunningJob( job ) );
}

RunningJob const &
get_local_running_job() {
	return * running_job_;
}

RunningJob::RunningJob( core::Size index ) : index_( index ) {}

core::Size
RunningJob::currently_running_job_index() const {
	return index_;
}

MTJob::MTJob() :
	done_( 0 ),
	index_( 0 ),
	rng_seed_( 0 )
{}

void
MTJob::go() {
	// 1. Initialize the RNG
	numeric::random::rg().set_seed( "mt19937", rng_seed_ );

	// 2. Initialize the thread-local running job
	RunningJob rj( index_ );
	set_running_job( rj );

	// 3. Get ready and run apply
	starttime_ = time(NULL);
	try {
		mover_->apply( *pose_ );
		mover_status_ = mover_->get_last_move_status();
	} catch (std::ios_base::failure & ex) {
		// OK: Here's the thing.  std::cout is a singular entity.  We can't go setting and unsetting the IO exit error code.

		std::cerr << "std::IO error detected... exiting..." << std::endl; // We can not longer use Tracer's at this point
		// Using pure exit instead of utility_exit_with_status to avoid infinite recursion
		std::exit( basic::options::option[basic::options::OptionKeys::out::std_IO_exit_error_code]());
	} catch (utility::excn::EXCN_BadInput& excn) {
		TR.Error
			<< "\n\n[ERROR] Exception caught by JobDistributor for job "
			<< output_name_ << excn
			<< std::endl;
		mover_status_ = protocols::moves::FAIL_BAD_INPUT;
	} catch (utility::excn::EXCN_Base& excn) {
		TR.Error
			<< "\n\n[ERROR] Exception caught by JobDistributor for job "
			<< output_name_ << excn
			<< std::endl;
		mover_status_ = protocols::moves::FAIL_DO_NOT_RETRY;
	}

	stoptime_ = time( NULL );

	// 4. inform the rest of the world that this thread has completed.
	++done_;
}

void MTJob::index( core::Size setting ) {
	index_ = setting;
}

void MTJob::rng_seed( int setting ) {
	rng_seed_ = setting;
}

void MTJob::mtjob_group( MTJobGroupOP setting ) {
	job_group_ = setting;
}

void MTJob::pose( core::pose::PoseOP pose ) {
	pose_ = pose;
}

void MTJob::mover( protocols::moves::MoverOP mover )
{
	mover_ = mover;
}

void MTJob::job( JobOP setting ) {
	job_ = setting;
}

void MTJob::output_name( std::string const & setting ) {
	output_name_ = setting;
}

core::pose::PoseOP
MTJob::pose()
{
	return pose_;
}

protocols::moves::MoverOP
MTJob::mover()
{
	return mover_;
}

protocols::moves::MoverStatus
MTJob::mover_status() const
{
	return mover_status_;
}

time_t
MTJob::starttime()
{
	return starttime_;
}

time_t
MTJob::stoptime()
{
	return stoptime_;
}


bool
MTJob::done() const
{
	return done_ == 1;
}

core::Size
MTJob::index() const
{
	return index_;
}

int
MTJob::rng_seed() const
{
	return rng_seed_;
}

MTJobGroupOP MTJob::mtjob_group() const
{
	return job_group_;
}

JobOP
MTJob::job() const
{
	return job_;
}

std::string const &
MTJob::output_name() const
{
	return output_name_;
}

core::Size & MTJob::n_retries()
{
	return job_group_->n_retries();
}

core::Size MTJob::n_retries() const
{
	return job_group_->n_retries();
}


MTJobGroup::MTJobGroup() :
	n_retries_( 0 )
{}

core::Size & MTJobGroup::n_retries()
{
	return n_retries_;
}

core::Size MTJobGroup::n_retries() const
{
	return n_retries_;
}


MultiThreadedJobDistributor::MultiThreadedJobDistributor() :
	JobDistributor(),
	max_n_running_threads_(0),
	sleep_min_( 10 ),
	sleep_max_( 2000 )
{
	max_n_running_threads_ = basic::options::option[ basic::options::OptionKeys::jd2::nthreads ];
}

MultiThreadedJobDistributor::~MultiThreadedJobDistributor()
{}

/// @details The main logic for the multithreaded job distributor. In the first loop
/// the MTJD tries to launch as many jobs as there are available threads to run until
/// they have all been launched, checking for jobs that have completed and writing
/// them out.  Once all jobs have been launched, then the MTJD waits until they have
/// all completed, writing out their results, before returning.
void MultiThreadedJobDistributor::go( protocols::moves::MoverOP mover )
{
	preliminaries();
	initialize_jobs_list();
	while ( jobs_remain() ) {
		if ( ! at_capacity() ) {
			launch_new_job( mover );
		} else {
			sleep_briefly();
		}
		check_for_completed_jobs();
	}

	while ( ! all_jobs_have_finished() ) {
		check_for_completed_jobs();
		sleep_briefly();
	}

}

void MultiThreadedJobDistributor::restart() {
	//	next_job_to_try_assigning_ = 1;
	mt_jobs_queue_.clear();
	//all_mt_jobs_.clear();
	mt_jobs_queue_.clear();
	running_jobs_.clear();
	JobDistributor::restart();
}

JobOP
MultiThreadedJobDistributor::current_job() const
{
	// With a one-job-per-thread model, we can respond to this request by figuring out
	// which job the current thread is running
	core::Size curr_job_index = get_local_running_job().currently_running_job_index();
	std::map< core::Size, MTJobOP >::const_iterator iter = running_jobs_.find( curr_job_index );
	assert( iter != running_jobs_.end() );
	return iter->second->job();
}

std::string
MultiThreadedJobDistributor::current_output_name() const
{
	core::Size curr_job_index = get_local_running_job().currently_running_job_index();
	std::map< core::Size, MTJobOP >::const_iterator iter = running_jobs_.find( curr_job_index );
	assert( iter != running_jobs_.end() );
	return iter->second->output_name();
}

core::Size
MultiThreadedJobDistributor::get_new_job_id()
{
	// This function should not be called.
	throw utility::excn::EXCN_Msg_Exception( "MultiThreadedJobDistributor should not have its get_new_job_id() function called" );
	return 0;
}

void
MultiThreadedJobDistributor::handle_interrupt()
{
	// ??!!
}


void
MultiThreadedJobDistributor::mark_current_job_id_for_repetition()
{
	curr_job_needs_requeuing_ = true;
}

void
MultiThreadedJobDistributor::remove_bad_inputs_from_job_list()
{
	curr_job_has_bad_inputs_ = true;
}


void MultiThreadedJobDistributor::preliminaries()
{
	// Permanently set the cout.e
	if (basic::options::option[basic::options::OptionKeys::out::std_IO_exit_error_code]()	> 0) {
		std::cout.exceptions(std::ios_base::badbit);
	}
}

int MultiThreadedJobDistributor::determine_starting_rng_seed()
{
	if ( basic::options::option[ basic::options::OptionKeys::run::constant_seed ] ) {
		return basic::options::option[ basic::options::OptionKeys::run::jran ];
	}
	return numeric::random::rg().get_seed();
}
void MultiThreadedJobDistributor::initialize_jobs_list() {
	JobsContainer const & parent_jobs_list( get_jobs() );
	//all_mt_jobs_.reserve( parent_job_list.size() );

	int seed = determine_starting_rng_seed();

	// Group together all the jobs with the same output_name() (i.e. that have different nstruct indices).
	// These jobs share a single counter for how many retries have been attempted for the job.
	MTJobGroupOP group( new MTJobGroup );
	mt_job_groups_.push_back( group );

	for ( core::Size ii = 1; ii <= parent_jobs_list.size(); ++ii ) {
		MTJobOP mtjob( new MTJob );
		mtjob->index( ii );
		mtjob->rng_seed( seed + ii );
		mtjob->job( parent_jobs_list[ ii ]->clone() );
		mtjob->output_name( job_outputter()->output_name( parent_jobs_list[ ii ] ));
		if ( ii != 1 && job_outputter()->output_name( parent_jobs_list[ ii-1 ] ) != job_outputter()->output_name( parent_jobs_list[ ii ] ) ) {
			group = MTJobGroupOP( new MTJobGroup );
			mt_job_groups_.push_back( group );
		}
		mtjob->mtjob_group( group );
		mt_jobs_queue_.push_back( mtjob );
		//all_mt_jobs_.push_back( mtjob );
	}

}

bool MultiThreadedJobDistributor::jobs_remain() const
{
	return ! mt_jobs_queue_.empty();
}

bool MultiThreadedJobDistributor::at_capacity() const
{
	assert( running_jobs_.size() <= max_n_running_threads_ );
	return running_jobs_.size() >= max_n_running_threads_;
}

void MultiThreadedJobDistributor::launch_new_job( protocols::moves::MoverOP mover )
{
	MTJobOP mtjob = mt_jobs_queue_.front();
	mt_jobs_queue_.pop_front();

	core::pose::PoseOP pose( new core::pose::Pose );
	job_inputter()->pose_from_job( *pose, mtjob->job() );

	if ( using_parser() ) {
		try {
			parser()->generate_mover_from_job( mtjob->job(), *pose, mover,
				mtjob->job()->nstruct_index() == 1,
				false /* do not store a modified Pose into the inner job */,
				true /* always generate a new mover */ );
		} catch ( utility::excn::EXCN_Base & excn ) {
			basic::Error()
				<< "ERROR: Exception caught by JobDistributor while trying to generate a Mover from job '"
				<< job_outputter()->output_name( mtjob->job() ) << "'" << std::endl
				<< excn
				<< std::endl;
			basic::Error()
				<< "Treating failure as bad input; canceling similar jobs"
				<< std::endl;
			remove_bad_inputs_from_job_list();
			job_failed(*pose, false);
		}
	} else {
		// uh oh! This is probably not going to work for you.
		// Soon, with serialization implemented, a deep copy of
		// the input mover can be created.
		mover = mover->clone();
	}

	// store the (possibly-modified-by-generate-mover-for-job) pose in the MTJob
	mtjob->pose( pose );
	mover->set_current_tag( job_outputter()->output_name( mtjob->job() ) );
	mtjob->mover( mover );

	running_jobs_[ mtjob->index() ] = mtjob;

	// Now launch a thread to run the job and detach from it, so that we
	// can exit this function while the other thread runs.

	std::thread newthread( [ mtjob ]() { mtjob->go(); } );
	newthread.detach();

}

void MultiThreadedJobDistributor::sleep_briefly() {
	if ( sleep_nsecs_ < sleep_max_ ) {
		sleep_nsecs_ *= 2;
	}

	std::chrono::nanoseconds dur( sleep_nsecs_ );
	std::this_thread::sleep_for( dur );
}

void MultiThreadedJobDistributor::reset_sleep_counters() {
	sleep_nsecs_ = sleep_min_;
}


void MultiThreadedJobDistributor::check_for_completed_jobs()
{
	// Iterate across the list of running jobs.  If any running job is reported as done,
	// then send its output to the job_outputter.
	for ( std::map< core::Size, MTJobOP >::iterator
			rj_iter     = running_jobs_.begin(),
			rj_iter_end = running_jobs_.end();
			rj_iter != rj_iter_end; /* no increment */ ) {
		std::map< core::Size, MTJobOP >::iterator rj_iter_next = rj_iter;
		++rj_iter_next;

		MTJobOP mtjob = rj_iter->second;
		if ( mtjob->done() ) {
			// We have found a completed job.  Rely on the parent class's write_output_from_job method.

			curr_job_needs_requeuing_ = false;
			curr_job_has_bad_inputs_ = false;
			JobDistributor::set_current_job_by_index( mtjob->index() );

			write_output_from_job(
				*mtjob->pose(), mtjob->mover(), mtjob->mover_status(),
				static_cast< core::Size > ( mtjob->stoptime() - mtjob->starttime() ),
				mtjob->n_retries() );

			// respond to any data-update requests made by the base class
			if ( curr_job_needs_requeuing_ ) {
				mt_jobs_queue_.push_front( mtjob );
			}
			//else {
			//	all_mt_jobs_[ mtjob->index() ] = 0; // Ensure this job will be deleted when mtjob goes out of scope.
			//}

			if ( curr_job_has_bad_inputs_ ) {
				// iterate across the queue and remove all jobs with the same output_tag:
				for ( std::list< MTJobOP >::iterator
						qiter     = mt_jobs_queue_.begin(),
						qiter_end = mt_jobs_queue_.end();
						qiter != qiter_end; /* no increment*/ ) {
					std::list< MTJobOP >::iterator qiter_next = qiter;
					++qiter_next;
					if ( (*qiter)->output_name() == mtjob->output_name() ) {
						mt_jobs_queue_.erase( qiter );
					}
					qiter = qiter_next;
				}
			}

			running_jobs_.erase( rj_iter );
			reset_sleep_counters();
		}

		rj_iter = rj_iter_next; // increment of iter occurs here.
	}
}

bool MultiThreadedJobDistributor::all_jobs_have_finished() const
{
	return running_jobs_.empty() && mt_jobs_queue_.empty();
}

}//jd2
}//protocols

#endif // CXX11
#endif // MULTI_THREADED
