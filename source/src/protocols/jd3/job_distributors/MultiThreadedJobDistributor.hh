// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MultiThreadedJobDistributor.hh
/// @brief  Job distributor that launches threads to carry out independant jobs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd2_MultiThreadedJobDistributor_hh
#define INCLUDED_protocols_jd2_MultiThreadedJobDistributor_hh

#ifdef MULTI_THREADED

// Unit headers
#include <protocols/jd2/MultiThreadedJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverStatus.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <map>
#include <string>

// Boost headers
#include <boost/detail/atomic_count.hpp>


namespace protocols {
namespace jd2 {

/// @brief Each worker thread will have a single instance of the RunningJob
/// that's used to track which job it is executing so that calls to
/// JobDistributor::current_job() can be appropriately handled.
class RunningJob : public utility::pointer::ReferenceCount
{
public:
	RunningJob( core::Size index );

	core::Size currently_running_job_index() const;

private:
	core::Size index_;
};


class MTJob : public utility::pointer::ReferenceCount
{
public:
	MTJob();

	void go();

	void index( core::Size setting );
	void rng_seed( int setting );
	void mtjob_group( MTJobGroupOP setting );
	void pose( core::pose::PoseOP setting );
	void mover( protocols::moves::MoverOP setting );
	void job( JobOP setting );
	void output_name( std::string const & setting );

	core::pose::PoseOP pose();
	protocols::moves::MoverOP mover();
	protocols::moves::MoverStatus mover_status() const;

	time_t starttime();
	time_t stoptime();
	bool done() const;

	core::Size index() const;
	int rng_seed() const;
	 MTJobGroupOP mtjob_group() const;

	JobOP job() const;
	std::string const & output_name() const;

	core::Size & n_retries();
	core::Size n_retries() const;


private:
	time_t starttime_;
	time_t stoptime_;
	boost::detail::atomic_count done_; // starts at 0, incremented to 1 when job is complete.
	core::Size index_;
	core::Size rng_seed_;
	core::pose::PoseOP pose_;
	protocols::moves::MoverOP mover_;
	protocols::moves::MoverStatus mover_status_;
	core::Size n_retries_;
	JobOP job_;
	MTJobGroupOP job_group_;
	std::string output_name_;
};

/// @brief A class holding data for a group of MTJobs that share the same output_name
/// i.e. that are different nstruct versions of the same InnerJob.
class MTJobGroup : public utility::pointer::ReferenceCount
{
public:
	MTJobGroup();
	core::Size & n_retries();
	core::Size n_retries() const;
private:
	core::Size n_retries_;
};

class MultiThreadedJobDistributor : public JobDistributor
{
protected:
	MultiThreadedJobDistributor();

public:
	///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
	///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
	virtual ~MultiThreadedJobDistributor();

	/// @brief The main workhorse of the JobDistibutor: it's handed a mover and then has to go
	/// and launch jobs that use that mover. For the MultiThreadedJobDistributor, it will need
	/// copies of that mover, so for now, the best way to use this JobDistributor is through
	/// RosettaScripts.
	virtual void go( protocols::moves::MoverOP mover );

	/// @brief Reset all internal data if the jobs are to be processed a second time (e.g. with a
	/// different mover).
	virtual void restart();

	/// @brief Movers may ask their controlling job distributor for information about the current job.
	/// They may also write information to this job for later output, though this use is now discouraged
	/// as the addition of the MultiplePoseMover now means that a single job may include several
	/// separate trajectories.
	virtual
	JobOP
	current_job() const;

	/// @brief Movers may ask their controlling job distributor for the output name as defined by the Job and JobOutputter.
	virtual
	std::string
	current_output_name() const;


	virtual core::Size get_new_job_id();
	virtual void handle_interrupt();


	friend class JobDistributorFactory; // calls protected ctor

protected:

	/// @brief The base class invokes this function during its write_output_for_job() whenever a
	/// job "soft-fails" and needs to be retried.  The MTJD handles this by re-queueing (at the
	/// beginning of the queue) the job that it's currently examining.
	virtual
	void
	mark_current_job_id_for_repetition();

	/// @brief The base class invokes this function during its write_output_for_job() whenever a
	/// job fails due to bad inputs. The MTJD responds to this call by removing from its job queue
	/// all jobs with the same output_name as the one that's currently being examined.
	virtual
	void
	remove_bad_inputs_from_job_list();


private:

	void preliminaries();

	int determine_starting_rng_seed();

	void initialize_jobs_list();

	/// @brief Have all jobs been launched?
	bool jobs_remain() const;

	/// @brief Are there as many running jobs as there is capacity?
	bool at_capacity() const;

	/// @brief Launch a new job using (a copy of) the input mover.
	void launch_new_job( protocols::moves::MoverOP mover );

	/// @brief Set the thread running the MTJD to sleep for a bit, increasing the
	/// amount of time it sleeps in each call, until reset_sleep_counters is called.
	void sleep_briefly();

	/// @brief Reset the amount of time that the MTJD spends sleeping before checking
	/// whether any of its running jobs have finished running.
	void reset_sleep_counters();

	/// @brief Look at the set of running jobs and if any of them have completed, then
	/// write their outputs, and mark the jobs as complete, if appropriate.  Also, reset
	/// the sleep counters.
	void check_for_completed_jobs();

	/// @brief Return true if there are no more jobs to launch and no more jobs still running.
	bool all_jobs_have_finished() const;

private:
	core::Size max_n_running_threads_;

	std::list< MTJobGroupOP >       mt_job_groups_;
	//utility::vector1< MTJobOP >     all_mt_jobs_;
	std::list< MTJobOP >            mt_jobs_queue_;
	std::map< core::Size, MTJobOP > running_jobs_;

	core::Size const sleep_min_;
	core::Size const sleep_max_;
	core::Size sleep_nsecs_;

	bool curr_job_needs_requeuing_;
	bool curr_job_has_bad_inputs_;

};

}//jd2
}//protocols

#endif // MULTI_THREADED

#endif //INCLUDED_protocols_jd2_MultiThreadedJobDistributor_HH
